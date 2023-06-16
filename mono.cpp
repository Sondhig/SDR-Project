/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "mode0.h"
#include "mode1.h"

int main(int argc, char *argv[]) {

    // VARIABLE INITIALIZATIONS
    int mode = 0;
    float rf_Fs_m0 = 2.4e6;
    float rf_Fs_m1 = 2.5e6;
    float rf_Fc = 100e3;
    float rf_taps = 151;
    float rf_decim = 10;
    float expand_factor = 24;
    float audio_Fs = 48e3;
    float mode1_Fs = 250e3 * expand_factor;
    float audio_decim_m0 = 5;
    float audio_decim_m1 = 125;
    float audio_taps_m0 = 151;
    float audio_taps_m1 = 151 * expand_factor;
    float audio_Fc = 16000;
    int block_size_m0 =
        1024 * rf_decim * audio_decim_m0 * 2; // defining block size
    int block_size_m1 = 1024 * rf_decim * audio_decim_m1 * 2;

    // If 1 argument is specified, set that value to mode
    if (argc == 2) {
        mode = atoi(argv[1]);
    }

    // If no arguments are specified or mode=0
    if (argc == 1 || mode == 0) {
        std::cerr << "Operating in mode 0" << std::endl;

        std::vector<float> rf_coeff;
        std::vector<float> audio_coeff;
        // computing the coefficients for convolution
        impulseResponseLPF(rf_Fs_m0, rf_Fc, rf_taps, rf_coeff);
        impulseResponseLPF(5 * audio_Fs, audio_Fc, audio_taps_m0, audio_coeff);
        std::vector<float> state_i_lpf_100k(rf_taps - 1, 0.0);
        std::vector<float> state_q_lpf_100k(rf_taps - 1, 0.0);
        std::vector<float> audio_zi(audio_taps_m0 - 1,
                                    0.0); // state used for mono block
        std::vector<float>
            audio_data; // Determine size instead of using push back
        std::vector<short int> audio_block_int;
        std::vector<float> prev_phase(2,
                                      0.0); // state used for our demod function
        std::vector<float> i_filt, q_filt;
        std::vector<float> i_data(block_size_m0 / 2, 0.0),
            q_data(block_size_m0 / 2, 0.0);
        int block_id = 0;

        std::vector<float> iq_data(block_size_m0);
        while (1) {
            // Reading data from stdin until end of stream
            readStdinBlockData(block_size_m0, block_id, iq_data);
            if ((std::cin.rdstate()) != 0) {
                break;
            }
            // Call mode0 function to process stdin data
            mode0(iq_data, block_size_m0, i_data, q_data, i_filt, q_filt,
                  rf_coeff, state_i_lpf_100k, state_q_lpf_100k, rf_decim,
                  prev_phase, audio_coeff, audio_zi, audio_decim_m0);

            block_id++;
            // Increment block
        }
    }

    else if (mode == 1) {
        std::cerr << "Operating in mode 1" << std::endl;

        std::vector<float> rf_coeff;
        std::vector<float> audio_coeff;
        // computing the coefficients for convolution
        impulseResponseLPF(rf_Fs_m1, rf_Fc, rf_taps, rf_coeff);
        impulseResponseLPF(mode1_Fs, audio_Fc, audio_taps_m1, audio_coeff);
        std::vector<float> state_i_lpf_100k(rf_taps - 1, 0.0);
        std::vector<float> state_q_lpf_100k(rf_taps - 1, 0.0);
        std::vector<float> audio_zi(audio_taps_m0 - 1,
                                    0.0); // state used for mono block
        std::vector<float>
            audio_data; // Determine size instead of using push back
        std::vector<short int> audio_block_int;
        std::vector<float> prev_phase(2,
                                      0.0); // state used for our demod function
        std::vector<float> i_filt, q_filt;
        std::vector<float> i_data(block_size_m1 / 2, 0.0),
            q_data(block_size_m1 / 2, 0.0);
        int block_id = 0;
        std::vector<float> iq_data(block_size_m1);

        while (1) {
            // Reading data from stdin until end of stream
            readStdinBlockData(block_size_m1, block_id, iq_data);
            if ((std::cin.rdstate()) != 0) {
                break;
            }
            // Call mode1 function to process stdin data
            mode1(iq_data, block_size_m1, i_data, q_data, i_filt, q_filt,
                  rf_coeff, state_i_lpf_100k, state_q_lpf_100k, rf_decim,
                  prev_phase, audio_coeff, audio_zi, audio_decim_m1,
                  expand_factor);
            block_id++;
            // Increment block
        }

    }

    else {
        // Argument entered is invalid
        std::cerr << "Invalid Mode" << std::endl;
        exit(1);
    }

    return 0;
}