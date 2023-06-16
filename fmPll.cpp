#include "dy4.h"
#include "iofunc.h"
#include <cmath>
#include <math.h>

void fmPll(std::vector<float> &pllIn, int freq, float Fs,
           std::vector<float> &state, float ncoScale, float phaseAdjust,
           float normBandwidth, std::vector<float> &ncoOut) {
    const float Cp = 2.666;
    const float Ci = 3.555;

    float Kp = normBandwidth * Cp;
    float Ki = normBandwidth * normBandwidth * Ci;

    ncoOut.resize(pllIn.size() + 1, 0.0);

    float integrator = state[0];
    float phaseEst = state[1];
    float feedbackI = state[2];
    float feedbackQ = state[3];
    ncoOut[0] = state[4];
    float trigOffset = state[5];

    for (int i = 0; i < pllIn.size(); i++) {

        // phase detector
        float errorI = pllIn[i] * feedbackI;
        float errorQ = pllIn[i] * -feedbackQ;

        float errorD = std::atan2(errorQ, errorI);

        // loop filter
        integrator = integrator + Ki * errorD;

        // update phase estimate
        phaseEst = phaseEst + Kp * errorD + integrator;

        // internal oscillator
        float trigArg = 2 * PI * (freq / Fs) * (trigOffset + i + 1) + phaseEst;
        feedbackI = std::cos(trigArg);
        feedbackQ = std::sin(trigArg);
        ncoOut[i + 1] = std::cos(trigArg * ncoScale + phaseAdjust);

        // new states
        state[0] = integrator;
        state[1] = phaseEst;
        state[2] = feedbackI;
        state[3] = feedbackQ;
        state[4] = ncoOut.back();
        state[5] = trigOffset + pllIn.size();
    }
}
