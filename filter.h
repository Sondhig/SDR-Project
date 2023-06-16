/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

/*a struct that contains the reference to the 1st element
of a vector and its size, so that we can pass by reference*/
struct blockVector {
    float *first_element;
    int size;
    blockVector(float *_first_element, int _size) {
        first_element = _first_element;
        size = _size;
    }
};

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseStereo(float, float, float, unsigned short int,
                           std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &,
                 const std::vector<float> &);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi);
void decimate(std::vector<float> &x, std::vector<float> &xDecimated, int n);
void Demod(const std::vector<float> &I, const std::vector<float> &Q,
           std::vector<float> &fm_demod, std::vector<float> &prev_phase);
void Demod(const std::vector<float> &I, const std::vector<float> &Q,
           std::vector<float> &fm_demod);
void arctanDemod(const std::vector<float> &I, const std::vector<float> &Q,
                 std::vector<float> &fm_demod, float prev_phase);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int n);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate, const int upsample_rate);
void Upsample(const std::vector<float> &x, std::vector<float> &xUpsampled,
              int n);

// These are the blockVector versions of the functions which we are overloading
void convolveFIR(std::vector<float> &y, blockVector x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate);
void convolveFIR(std::vector<float> &y, blockVector x,
                 const std::vector<float> &h, std::vector<float> &zi,
                 const int downsample_rate, const int upsample_rate);
void Upsample(blockVector x, std::vector<float> &xUpsampled, int n);

#endif // DY4_FILTER_H
