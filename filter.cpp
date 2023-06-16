/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "filter.h"

#include <math.h>

#include <cmath>

#include "dy4.h"
#include "iofunc.h"
// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps,
                        std::vector<float> &h) {
    // allocate memory for the impulse response
    h.resize(num_taps, 0.0);
    float norm = 2 * Fc / Fs;            // initialize norm
    for (int i = 0; i < num_taps; i++) { // iterate over number of taps
        if (i == ((num_taps - 1) / 2)) { // calculate value of h[i] using the
                                         // impulse response summation equation
            h[i] = norm;
        } else {
            h[i] = norm * ((std::sin(PI * norm * (i - (num_taps - 1) / 2))) /
                           (PI * norm * (i - (num_taps - 1) / 2)));
        }
        h[i] *= pow((std::sin(i * PI / num_taps)), 2);
    }
}
void impulseResponseStereo(float Fb, float Fe, float Fs,
                           unsigned short int num_taps, std::vector<float> &h) {
    // allocate memory for the impulse response
    h.resize(num_taps, 0.0);
    float norm_cent = ((Fe + Fb) / 2) / (Fs / 2); // initialize norm
    float norm_pass = (Fe - Fb) / (Fs / 2);
    for (int i = 0; i < num_taps; i++) { // iterate over number of taps
        if (i == ((num_taps - 1) / 2)) { // calculate value of h[i] using the
                                         // impulse response summation equation
            h[i] = norm_pass;
        } else {
            h[i] =
                norm_pass *
                ((std::sin(PI * (norm_pass / 2) * (i - (num_taps - 1) / 2))) /
                 (PI * (norm_pass / 2) * (i - (num_taps - 1) / 2)));
        }
        h[i] *= std::cos(i * PI * norm_cent);
        h[i] *= pow((std::sin(i * PI / num_taps)), 2);
    }
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h) {
    // allocate memory for the output (filtered) data
    y.resize(x.size() + h.size() - 1, 0.0);
    for (int i = 0; i < x.size(); i++) { // loop through the size of x
        float sum = 0;
        for (int j = 0; j < h.size(); j++) { // then loop through the size of h
            if ((i - j) >= 0) { // only convolute if the indexes are both
                                // positive (otherwise for block processing
                                // would need previous block data)
                sum += h[j] * x[i - j];
            }
        }
        y[i] = sum; // add value to y index
    }
}

// overloaded convolveFIR function with zi as a parameter
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi) {
    // allocate memory for the output (filtered) data
    y.resize(x.size(), 0.0);
    // std::cerr << "x.size: " << x.size() << " \n"; //151
    // std::cerr << "zi.size: " << zi.size() << " \n"; //150
    // std::cerr << "h.size: " << h.size() << " \n"; //51200
    for (int i = 0; i < x.size(); i++) { // loop through the size of x
        float sum = 0;
        for (int j = 0; j < h.size(); j++) { // then loop through the size of h
            // std::cerr << "i-j: " << i-j << std::endl;
            if ((i - j) >= 0) { // if the indexes are both positive
                sum += h[j] * x[i - j];
            } else { // if not then use values from zi
                // std::cerr << "index: " << zi.size()+(i-j) << "\n";
                sum += h[j] * zi[zi.size() + (i - j)];
            }
        }
        y[i] = sum; // add value to y index
    }
    // calculating the new zi from x
    auto first = x.cbegin() + (x.size() - h.size() + 1);
    auto last = x.cbegin() + (x.size());
    std::vector<float> newZi(first, last);
    zi = newZi;
}

// faster convolveFIR function with the decimation that follows it
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate) {
    // allocate memory for the output (filtered) data
    int k = 0;
    y.resize(x.size() / downsample_rate, 0.0);
    for (int i = 0; i < x.size();
         i += downsample_rate) { // loop through the size of x
        // if (i%n==0) { //if it will be kept at decimation
        /*since block_size is the multiple of rf_decim and audio_decim,
         i%n==0 can be used. Otherwise, (i+block_size*block_count)%n==0
         should be used*/
        float sum = 0;
        for (int j = 0; j < h.size(); j++) { // then loop through the size of h
            // std::cerr << "i-j: " << i-j << std::endl;
            if ((i - j) >= 0) { // if the indexes are both positive
                sum += h[j] * x[i - j];
            } else { // if not then use values from zi
                // std::cerr << "index: " << zi.size()+(i-j) << "\n";
                sum += h[j] * zi[zi.size() + (i - j)];
            }
        }
        y[k] = sum; // add value to y index
        k++;
        // }
    }
    // calculating the new zi from x
    auto first = x.cbegin() + (x.size() - h.size() + 1);
    auto last = x.cbegin() + (x.size());
    std::vector<float> newZi(first, last);
    zi = newZi;
}

void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate, const int upsample_rate) {
    // n is the decimation factor
    // allocate memory for the output (filtered) data
    int k = 0;
    // y.resize(x.size()/downsample_rate, 0.0);
    for (int i = 0; i < x.size();
         i += downsample_rate) { // loop through the size of x
        float sum = 0;
        int start = (k * downsample_rate) % upsample_rate;
        for (int j = start; j < h.size();
             j += upsample_rate) { // then loop through the size of h
            if ((i - j) >= 0) {    // if the indexes are both positive
                sum += h[j] * x[i - j];
            } else { // if not then use values from zi
                sum += h[j] * zi[zi.size() + (i - j)];
            }
        }
        y[k] = sum * upsample_rate; // add value*upsample rate to y index
        k++;
    }

    // calculating the new zi from x
    auto first = x.cbegin() + (x.size() - h.size() + 1);
    auto last = x.cbegin() + (x.size());
    std::vector<float> newZi(first, last);
    zi = newZi;
}

void decimate(std::vector<float> &x, std::vector<float> &xDecimated, int n) {
    float length = x.size();
    int newLength = ceil(length / n);
    // std::cerr << "New Length " << newLength << "\n";

    xDecimated.resize(newLength, static_cast<float>(0));

    int j = 0;
    for (int i = 0; i < length; i += n) {
        xDecimated[j++] = x[i];
    }
}
// without prev_phase parameter
void Demod(const std::vector<float> &I, const std::vector<float> &Q,
           std::vector<float> &fm_demod) {
    fm_demod.resize(I.size(), 0.0);
    printRealVector(fm_demod);
    for (auto k = 0; k < I.size(); k++) {
        if (k == 0) {
            fm_demod[k] = 0;
        } else {
            fm_demod[k] =
                (I[k] * (Q[k] - Q[k - 1]) - Q[k] * (I[k] - I[k - 1])) /
                (pow(I[k], 2) + pow(Q[k], 2));
        }
    }
}

void Demod(const std::vector<float> &I, const std::vector<float> &Q,
           std::vector<float> &fm_demod, std::vector<float> &prev_phase) {
    fm_demod.resize(I.size(), 0.0);
    for (auto k = 0; k < I.size(); k++) {
        auto denominator = (pow(I[k], 2) + pow(Q[k], 2));
        if (k == 0) {
            if (denominator == 0) {
                fm_demod[k] = 0;
            } else {
                fm_demod[k] = (I[k] * (Q[k] - prev_phase[1]) -
                               Q[k] * (I[k] - prev_phase[0])) /
                              denominator;
            }
        } else {
            fm_demod[k] =
                (I[k] * (Q[k] - Q[k - 1]) - Q[k] * (I[k] - I[k - 1])) /
                denominator;
        }
    }
    prev_phase[0] = I[I.size() - 1];
    prev_phase[1] = Q[Q.size() - 1];
}

// Implementing included python arctan demod
void arctanDemod(const std::vector<float> &I, const std::vector<float> &Q,
                 std::vector<float> &fm_demod, float prev_phase) {
    fm_demod.resize(I.size(), 0.0);
    float current_phase, diff;
    for (auto k = 0; k < I.size(); k++) {
        current_phase = atan2(Q[k], I[k]);
        diff = current_phase - prev_phase;
        // if (k<10) {
        // std::cerr << "diff for "<<k<< "is "<<diff<<std::endl;}
        // the two if statements below are like np.unwrap
        // confines the diff within -pi and pi rad
        while (diff > PI) {
            diff = diff - 2 * PI;
        }
        while (diff < -PI) {
            diff = diff + 2 * PI;
        }
        // if (k<15) {
        // std::cerr << "demod for "<<k<< "is "<<diff<<std::endl;
        // }
        fm_demod[k] = diff;
        prev_phase = current_phase;
    }
}

// void Upsample(std::vector<float> &x,std::vector<float> &xUpsampled, int
// n){//add n-1 zeros 	float length = x.size(); 	int newLength =
// ceil(length * n); 	xUpsampled.resize(newLength, static_cast<float>(0));
// int
// j=0; 	for(int i=0;i<length;i++){ 		xUpsampled[j++]=x[i];
// for(int k=0;k<n-1;k++){ 			xUpsampled[j++]=0;
// 		}
// 	}
// }

// New Upsample with less operations

void Upsample(const std::vector<float> &x, std::vector<float> &xUpsampled,
              int n) { // add n-1 zeros
    // Requires a vector of proper size to be passed in
    float length = x.size();
    int newLength = length * n;
    int j = 0;
    for (int i = 0; i < newLength; i += n) {
        xUpsampled[i] = x[j++];
    }
}

// Start of the overloaded fxns for blockVectors
void convolveFIR(std::vector<float> &y, blockVector x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate, const int upsample_rate) {
    // n is the decimation factor
    // allocate memory for the output (filtered) data
    int k = 0;
    // y.resize(x.size()/downsample_rate, 0.0);
    for (int i = 0; i < x.size;
         i += downsample_rate) { // loop through the size of x
        float sum = 0;
        int start = (k * downsample_rate) % upsample_rate;
        for (int j = start; j < h.size();
             j += upsample_rate) { // then loop through the size of h
            if ((i - j) >= 0) {    // if the indexes are both positive
                sum += h[j] * (*(x.first_element + i - j));
            } else { // if not then use values from zi
                sum += h[j] * zi[zi.size() + (i - j)];
            }
        }
        y[k] = sum * upsample_rate; // add value*upsample rate to y index
        k++;
    }

    // calculating the new zi from x
    auto first = x.first_element + (x.size - h.size() + 1);
    auto last = x.first_element + (x.size);
    std::vector<float> newZi(first, last);
    zi = newZi;
}

// This is the one using our new structure
void convolveFIR(std::vector<float> &y, blockVector x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate) {
    // n is the decimation factor
    // allocate memory for the output (filtered) data
    int k = 0;
    y.resize(x.size / downsample_rate, 0.0);
    for (int i = 0; i < x.size;
         i += downsample_rate) { // loop through the size of x
        // if (i%n==0) { //if it will be kept at decimation
        /*since block_size is the multiple of rf_decim and audio_decim,
         i%n==0 can be used. Otherwise, (i+block_size*block_count)%n==0
         should be used*/
        float sum = 0;
        for (int j = 0; j < h.size(); j++) { // then loop through the size of h
            // std::cerr << "i-j: " << i-j << std::endl;
            if ((i - j) >= 0) { // if the indexes are both positive
                sum += h[j] * (*(x.first_element + i - j));
            } else { // if not then use values from zi
                // std::cerr << "index: " << zi.size()+(i-j) << "\n";
                sum += h[j] * zi[zi.size() + (i - j)];
            }
        }
        y[k] = sum; // add value to y index
        k++;
        // }
    }
    // calculating the new zi from x
    auto first = x.first_element + (x.size - h.size() + 1);
    auto last = x.first_element + (x.size);
    std::vector<float> newZi(first, last);
    zi = newZi;
}

// This is for our new data structure
void Upsample(blockVector x, std::vector<float> &xUpsampled,
              int n) { // add n-1 zeros
    // Requires a vector of proper size to be passed in
    float length = x.size;
    int newLength = length * n;
    int j = 0;
    for (int i = 0; i < newLength; i += n) {
        xUpsampled[i] = *(x.first_element + j++);
    }
}
