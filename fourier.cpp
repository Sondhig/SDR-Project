/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "fourier.h"
#include "dy4.h"
#include <math.h>

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
    Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
    for (auto m = 0; m < Xf.size(); m++) {
        for (auto k = 0; k < x.size(); k++) {
            std::complex<float> expval(0, -2 * PI * (k * m) / x.size());
            Xf[m] += x[k] * std::exp(expval);
        }
    }
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf,
                            std::vector<float> &Xmag) {
    // only the positive frequencies
    Xmag.resize(Xf.size(), static_cast<float>(0));
    for (auto i = 0; i < Xf.size(); i++) {
        Xmag[i] = std::abs(Xf[i]) / Xf.size();
    }
}

// add your own code to estimate the PSD
void estimatePSD(const std::vector<float> &samples, const float Fs,
                 std::vector<float> &freq, std::vector<float> &psd_est) {
    int freq_bins = NFFT;
    float df = Fs / freq_bins;

    for (float value = 0; value < (Fs / 2); value += df)
        freq.push_back(value);
    std::vector<float> hann;
    hann.resize(freq_bins, 0.0);
    for (int i = 0; i < hann.size(); i++) {
        hann[i] = pow(sin(i * PI / freq_bins), 2);
    }
    std::vector<float> psd_list;
    int no_segments = floor(samples.size() / freq_bins);

    for (int i = 0; i < no_segments; i++) {
        std::vector<float>
            windowed_samples; // initializing windowed_samples array
        windowed_samples.resize(freq_bins, 0.0);
        for (int j = 0; j < freq_bins; j++) {
            windowed_samples[j] =
                samples[i * freq_bins + j] *
                hann[j]; // pointwise multiplication, applying Hann window
        }
        std::vector<std::complex<float>> Xf;
        DFT(windowed_samples, Xf);
        std::vector<std::complex<float>> Xf_split(
            Xf.cbegin(),
            Xf.cbegin() +
                int(freq_bins / 2)); // keeping only the positive half freqs
        float psd_seg;
        for (int j = 0; j < Xf_split.size(); j++) {
            psd_seg = std::pow(std::abs(Xf_split[j]),
                               2); // abs(Xf_split)^2,=Re^2+Im^2
            psd_seg *= 2 / (Fs * freq_bins / 2);
            psd_seg = 10 * std::log10(psd_seg);
            psd_list.push_back(psd_seg);
        }
    }
    psd_est.resize(int(freq_bins / 2), 0.0);
    for (int k = 0; k < int(freq_bins / 2); k++) {
        for (int l = 0; l < no_segments; l++) {
            psd_est[k] += psd_list[k + l * int(freq_bins / 2)];
        }
        psd_est[k] /= no_segments;
    }
}