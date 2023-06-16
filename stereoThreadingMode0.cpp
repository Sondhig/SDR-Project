#include <chrono>
#include <condition_variable>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>

#include "dy4.h"
#include "filter.h"
#include "fmPll.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <unistd.h>

// VARIABLE INITIALIZATIONS

#define QUEUE_BLOCKS 10
#define DEMOD_SIZE 5120 // size of one block after demodulation

float rf_Fs = 2.4e6;
float rf_Fc = 100e3;
float rf_taps = 151;
float rf_decim = 10;

float audio_Fs = 48e3;
float audio_decim = 5;
float audio_taps = 151;
float audio_Fc = 16000;

float pll_bpf_fb = 18500;
float pll_bpf_fe = 19500;
float stereo_bpf_fb = 22000;
float stereo_bpf_fe = 54000;

std::vector<float> rf_coeff;
std::vector<float> audio_coeff;
std::vector<float> pilot_tone_coeff;
std::vector<float> stereo_coeff;

std::vector<float> state_i_lpf_100k(rf_taps - 1, 0.0);
std::vector<float> state_q_lpf_100k(rf_taps - 1, 0.0);

std::vector<float> audio_zi(audio_taps - 1, 0.0); // state used for mono block
// states for stereo
std::vector<float> pilot_tone_zi(audio_taps - 1, 0.0);
std::vector<float> stereo_zi(audio_taps - 1, 0.0);
std::vector<float> mixer_zi(audio_taps - 1, 0.0);
std::vector<float> state_save{0.0, 0.0, 1.0, 0.0, 1.0, 0.0}; // used for pll

std::vector<float> audio_data;
std::vector<short int> audio_block_int;

std::vector<short int> stereo_data;
std::vector<float> prev_phase(2, 0.0);
std::vector<float> i_filt, q_filt, pilot_tone_filt, stereo_filt, mixer_block,
    mixer, pll;

std::queue<blockVector>
    demod_queue; // Queue of custom data type defined in filter.h

// defining the queue to hold demodulated blocks
static std::vector<std::vector<float>>
    queue_block(QUEUE_BLOCKS, std::vector<float>(DEMOD_SIZE, 0.0));

// rf thread. input: 2.4MHz 8bit unsigned int; output: 240kHz 32bit signed float
int rfThread(std::vector<float> rf_coeff, int block_size,
             std::vector<float> i_data, std::vector<float> q_data,
             std::mutex &qtex, std::condition_variable &cvar, int block_id) {

    int queue_entry = block_id % QUEUE_BLOCKS;

    // convolution and decimation. cut-off freq=100k, decimation factor 10
    convolveFIR(i_filt, i_data, rf_coeff, state_i_lpf_100k, rf_decim);
    convolveFIR(q_filt, q_data, rf_coeff, state_q_lpf_100k, rf_decim);

    std::vector<float> audio_filt;
    // for stereo
    std::vector<float> stereo_block;

    // demodulation using the demod function we wrote
    Demod(i_filt, q_filt, queue_block[queue_entry], prev_phase);

    // locking mutex
    std::unique_lock<std::mutex> my_lock(qtex);

    if (demod_queue.size() == QUEUE_BLOCKS) { // if queue full, wait
        cvar.wait(my_lock);
    }

    // blockVector constructor takes a pointer to the 0th element of the vector,
    // and demod size
    blockVector fm_demod((float *)&queue_block[queue_entry][0], DEMOD_SIZE);

    // Since this is a object we can just push it on the queue.
    demod_queue.push(fm_demod);
    // PUSH fm_demod TO QUEUE

    // unlocking mutex
    my_lock.unlock();
    cvar.notify_one();
    return 1;
}

// audio thread. input:240kHz 32bit signed float, output: 16bit signed int
// stereo
int audioThread(std::mutex &qtex, std::condition_variable &cvar) {
    std::vector<float> audio_block;

    // locking mutex
    std::unique_lock<std::mutex> my_lock(qtex);
    if (demod_queue.empty()) { // if queue empty, wait
        cvar.wait(my_lock);
    }

    // Since it is an object it can come right off the queue.
    blockVector fm_demod = demod_queue.front();

    demod_queue.pop();

    // unlocking mutex
    my_lock.unlock();
    cvar.notify_one();

    // convoluting to extract pilot tone and stereo data by applying the BPFs
    convolveFIR(pilot_tone_filt, fm_demod, pilot_tone_coeff, pilot_tone_zi,
                1.0);
    convolveFIR(stereo_filt, fm_demod, stereo_coeff, stereo_zi, 1.0);

    // PLL
    fmPll(pilot_tone_filt, 19000, 5 * audio_Fs, state_save, 2.0, 0.0, 0.01,
          pll);

    // mixer
    mixer.resize(stereo_filt.size(), 0.0);
    for (int i = 0; i < stereo_filt.size(); i++) {
        mixer[i] = 2 * pll[i] * stereo_filt[i];
        // 2* is because the demodulated signal only has 1/2 original
        // strength
    }

    // stereo: convolution Fc=16k with decimation by 5, just like mono
    convolveFIR(mixer_block, mixer, audio_coeff, mixer_zi, audio_decim);

    // mono: convolution. cut-off freq=16k with decimation by 5
    convolveFIR(audio_block, fm_demod, audio_coeff, audio_zi, audio_decim);

    // stereo 2-d array with left and right channels
    // initializing 16 bit signed int left and right channels
    std::vector<float> audio_block_left(audio_block.size());
    std::vector<float> audio_block_right(audio_block.size());
    // initializing stereo as 2x size of left/right
    std::vector<short int> stereo_block_int(audio_block.size() * 2);

    // getting the final stereo data and converting to 16bit signed int
    for (int i = 0; i < audio_block.size(); i++) {
        // stereo combining
        audio_block_left[i] = (audio_block[i] + mixer_block[i]) / 2;
        audio_block_right[i] = (audio_block[i] - mixer_block[i]) / 2;
        // if-else statements are for typecasting to short int
        if (std::isnan(audio_block_left[i])) {
            stereo_block_int[i * 2] = 0;
        } else {
            stereo_block_int[i * 2] =
                static_cast<short int>(audio_block_left[i] * 16384);
        }
        if (std::isnan(audio_block_right[i])) {
            stereo_block_int[i * 2 + 1] = 0;
        } else {
            stereo_block_int[i * 2 + 1] =
                static_cast<short int>(audio_block_right[i] * 16384);
        }
    }

    // writing to stdout to be heard by aplay
    fwrite(&stereo_block_int[0], sizeof(short int), stereo_block_int.size(),
           stdout);
    return 1;
}

int main() {
    // defining parameters
    // computing the coefficients for convolution
    impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
    impulseResponseLPF(5 * audio_Fs, audio_Fc, audio_taps, audio_coeff);
    // stereo coefficients for BPF convolution
    impulseResponseStereo(pll_bpf_fb, pll_bpf_fe, 5 * audio_Fs, audio_taps,
                          pilot_tone_coeff);
    impulseResponseStereo(stereo_bpf_fb, stereo_bpf_fe, 5 * audio_Fs,
                          audio_taps, stereo_coeff);

    int block_size = 1024 * rf_decim * audio_decim * 2; // defining block size

    // defining mutex and conditional variable
    std::mutex qtex;
    std::condition_variable cvar;

    // declaring some variables that will be used inside the while loop

    std::vector<float> i_data(block_size / 2, 0.0), q_data(block_size / 2, 0.0);

    int block_id = 0; // counts how many blocks have been processed
    std::vector<float> iq_data(block_size);

    while (1) {
        std::vector<float> iq_data(block_size);
        // Reading data from stdin until end of stream
        readStdinBlockData(block_size, block_id, iq_data);
        if ((std::cin.rdstate()) != 0) {
            break;
        }

        // Separating data from iq data to i and q components
        for (int i = 0; i < block_size / 2; i++) {
            i_data[i] = iq_data[2 * i];
            q_data[i] = iq_data[2 * i + 1];
        }

        // ---------- THREADS BEGIN -----------------------
        std::thread t_rf = std::thread(
            rfThread, std::ref(rf_coeff), std::ref(block_size),
            std::ref(i_data), std::ref(q_data), std::ref(qtex), std::ref(cvar),
            std::ref(block_id)); // spawn new thread that calls rfThread() and
                                 // passing in rf_coeff, block_size, i_data,
                                 // q_data, qtex, cvar, and block_id

        std::thread t_audio = std::thread(
            audioThread, std::ref(qtex),
            std::ref(cvar)); // spawn new thread that calls audioThread() and
                             // passing in qtex, cvar

        t_rf.join();
        t_audio.join();

        // ---------- THREADS END -----------------------

        block_id++;
    }
    return 0;
}