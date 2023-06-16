#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

void mode0(const std::vector<float> &iq_data, int block_size,
           std::vector<float> &i_data, std::vector<float> &q_data,
           std::vector<float> &i_filt, std::vector<float> &q_filt,
           std::vector<float> &rf_coeff, std::vector<float> &state_i_lpf_100k,
           std::vector<float> &state_q_lpf_100k, float rf_decim,
           std::vector<float> &prev_phase, std::vector<float> &audio_coeff,
           std::vector<float> &audio_zi, float audio_decim) {

    // Seperating iq_data into i and q vectors
    for (int i = 0; i < block_size / 2; i++) {
        i_data[i] = iq_data[2 * i];
        q_data[i] = iq_data[2 * i + 1];
    }

    // convolution. cut-off freq=100k, decimating by rf_decim
    convolveFIR(i_filt, i_data, rf_coeff, state_i_lpf_100k, rf_decim);
    convolveFIR(q_filt, q_data, rf_coeff, state_q_lpf_100k, rf_decim);
    std::vector<float> fm_demod, audio_filt, audio_block;

    // demodulation using the demod function we wrote
    Demod(i_filt, q_filt, fm_demod, prev_phase);

    // ------ABOVE: RF FRONTEND - BELOW: MONO BLOCK------

    convolveFIR(audio_block, fm_demod, audio_coeff, audio_zi,
                audio_decim); // convolution. cut-off freq=16k

    std::vector<short int> audio_data_int(
        audio_block.size()); // short int = 16 bit int

    for (int i = 0; i < audio_block.size(); i++) {
        if (std::isnan(audio_block[i]))
            audio_data_int[i] = 0;
        else
            audio_data_int[i] = static_cast<short int>(
                audio_block[i] * 16384); // Map to a 16 bit signed int
    }
    // Write to standard out
    fwrite(&audio_data_int[0], sizeof(short int), audio_data_int.size(),
           stdout);
}