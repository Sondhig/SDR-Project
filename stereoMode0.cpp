#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>

#include "dy4.h"
#include "filter.h"
#include "fmPll.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
// using namespace std;

int main() {
    // defining parameters
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

    // computing the coefficients for convolution
    impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
    impulseResponseLPF(5 * audio_Fs, audio_Fc, audio_taps, audio_coeff);
    // stereo coefficients for BPF convolution
    impulseResponseStereo(pll_bpf_fb, pll_bpf_fe, 5 * audio_Fs, audio_taps,
                          pilot_tone_coeff);
    impulseResponseStereo(stereo_bpf_fb, stereo_bpf_fe, 5 * audio_Fs,
                          audio_taps, stereo_coeff);

    int block_size = 1024 * rf_decim * audio_decim * 2; // defining block size

    // initializing the states (i.e.,zi parameters) for I and Q
    std::vector<float> state_i_lpf_100k;
    std::vector<float> state_q_lpf_100k;
    state_i_lpf_100k.resize(rf_taps - 1, 0.0);
    state_q_lpf_100k.resize(rf_taps - 1, 0.0);

    std::vector<float> audio_zi; // state used for mono block
    audio_zi.resize(audio_taps - 1, 0.0);
    // states for stereo
    std::vector<float> pilot_tone_zi;
    std::vector<float> stereo_zi;
    std::vector<float> mixer_zi;
    std::vector<float> state_save; // used for pll
    pilot_tone_zi.resize(audio_taps - 1, 0.0);
    stereo_zi.resize(audio_taps - 1, 0.0);
    mixer_zi.resize(audio_taps - 1, 0.0);
    state_save = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0}; // not sure if this is
                                                 // correct!

    std::vector<float> audio_data; // Determine size instead of using push back
    std::vector<short int> audio_block_int;

    // final output vector for stereo
    std::vector<short int> stereo_data;

    // COMMENT OUT one of the next two lines based on the demod function used
    std::vector<float> prev_phase(2, 0.0); // state used for our demod function
    // float prev_phase=0.0; //state used for arctanDemod from python

    // for measuring run time
    auto t_start = std::chrono::high_resolution_clock::now();

    // declaring some variables that will be used inside the while loop
    std::vector<float> i_filt, q_filt, pilot_tone_filt, stereo_filt,
        mixer_block;
    std::vector<float> i_data, q_data, mixer, pll;
    i_data.resize(block_size / 2, 0.0);
    q_data.resize(block_size / 2, 0.0);
    int block_id = 0; // counts how many blocks have been processed

    while (1) {
        std::cerr << "starting block " << block_id << std::endl;
        std::vector<float> iq_data(block_size);
        readStdinBlockData(block_size, block_id, iq_data);
        if ((std::cin.rdstate()) != 0) {
            // std::cerr << "should be exiting" << std::endl;
            break;
        }
        for (int i = 0; i < block_size / 2; i++) {
            i_data[i] = iq_data[2 * i];
            q_data[i] = iq_data[2 * i + 1];
        }
        // convolution. cut-off freq=100k
        // convolveFIRwDecim is faster than convolveFIR by taking decimation
        // into consideration
        convolveFIR(i_filt, i_data, rf_coeff, state_i_lpf_100k, rf_decim);
        convolveFIR(q_filt, q_data, rf_coeff, state_q_lpf_100k, rf_decim);

        std::vector<float> fm_demod, audio_filt, audio_block;
        // for stereo
        std::vector<float> stereo_block;

        // demodulation
        // using the demod function we wrote
        Demod(i_filt, q_filt, fm_demod, prev_phase);

        // ------BELOW: STEREO BLOCK; ABOVE: RF FRONTEND------

        // applying stereo band pass filters
        convolveFIR(pilot_tone_filt, fm_demod, pilot_tone_coeff, pilot_tone_zi);
        convolveFIR(stereo_filt, fm_demod, stereo_coeff, stereo_zi);

        // PLL, Should be correct
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

        if (block_id == 3) {
            std::vector<float> freq0, freq1, freq2;
            std::vector<float> psd_est_demod, psd_est_filt, psd_est_block;
            estimatePSD(fm_demod, (rf_Fs / rf_decim) / 1e3, freq0,
                        psd_est_demod);
            logVector("psd_est_demod", freq0, psd_est_demod);

            estimatePSD(pll, (rf_Fs / rf_decim) / 1e3, freq1, psd_est_filt);
            logVector("psd_est_filt", freq1,
                      psd_est_filt); // dont change file names

            estimatePSD(pilot_tone_filt, (rf_Fs / rf_decim) / 1e3, freq2,
                        psd_est_block);
            logVector("psd_est_block", freq2, psd_est_block);
        }

        // stereo 2-d array with left and right channels
        // initializing 16 bit signed int left and right channels
        // std::vector<short int> audio_block_left_int(audio_block.size());
        // std::vector<short int> audio_block_right_int(audio_block.size());
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
        // NOT SURE IF THE WAY BELOW OF WRITING INTO STEREO IS CORRECT
        // but looking from wavio.py this interleaving way is how we do it
        // stereo_block_int[i*2]=audio_block_left_int[i];
        // stereo_block_int[i*2+1]=audio_block_right_int[i];
        // stereo_data.insert(stereo_data.end(), stereo_block_int.begin(),
        // stereo_block_int.end());
        fwrite(&stereo_block_int[0], sizeof(short int), stereo_block_int.size(),
               stdout);
        block_id++;
    }

    // measuring run time
    auto t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> run_time = t_stop - t_start;
    std::cerr << "time taken in milliseconds: " << run_time.count();

    writeBinData("../data/testoutput.bin", stereo_data);
    std::cerr << "Finished block processing" << std::endl;

    return 0;
}