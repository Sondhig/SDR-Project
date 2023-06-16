/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// add headers as needed
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// declaration of a function prototypes
void printRealVector(const std::vector<float> &);

void printComplexVector(const std::vector<std::complex<float>> &);

void readBinData(const std::string, std::vector<float> &);

void writeBinData(const std::string, const std::vector<short int> &);

void ourWriteBinData(const std::string out_fname, const std::vector<float> &x);

void readStdinBlockData(int num_samples, unsigned int block_id,
                        std::vector<float> &block_data);

#endif // DY4_IOFUNC_H
