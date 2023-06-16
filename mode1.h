#include <vector>

void mode1(const std::vector<float> &iq_data, int block_size,
           std::vector<float> &i_data, std::vector<float> &q_data,
           std::vector<float> &i_filt, std::vector<float> &q_filt,
           std::vector<float> &rf_coeff, std::vector<float> &state_i_lpf_100k,
           std::vector<float> &state_q_lpf_100k, float rf_decim,
           std::vector<float> &prev_phase, std::vector<float> &audio_coeff,
           std::vector<float> &audio_zi, float audio_decim,
           float expand_factor);