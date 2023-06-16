#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan
from fmSupportLib import Demod

np.set_printoptions(threshold=np.inf)
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
audio_taps = 151
audio_Fc = 16000
# add other settings for audio, like filter taps, ...


def our_fir(fc, fs, Ntaps):
    # FIR=finite impulse response
    h = np.zeros(Ntaps)
    Norm = 2 * fc / fs
    for i in range(Ntaps):
        if i == (Ntaps - 1) / 2:
            h[i] = Norm
        else:
            h[i] = Norm * (
                (math.sin(math.pi * Norm * (i - (Ntaps - 1) / 2)))
                / (math.pi * Norm * (i - (Ntaps - 1) / 2))
            )
        h[i] = h[i] * math.pow((math.sin(i * math.pi / Ntaps)), 2)
    return h


def our_lfilter(b, a, x, N_taps, zi=None):
    y = np.zeros(len(x))  # initializing output list
    for i in range(len(x)):  # iterate for each y index
        sum = 0
        for j in range(N_taps):  # do the multiplication for convolution N_taps times
            if i - j >= 0:  # if it is calculated within the current block
                sum += b[j] * x[i - j]  # take the values from the input x
            elif not (zi is None):
                sum += b[j] * zi[i - j]  # otherwise take it from the zi list
        y[i] = sum  # save the value of sum the the index of the output Y
    if zi is None:
        return y
    return (
        y,
        x[len(x) - N_taps + 1 :],
    )  # return output Y and the last Ntaps-1 values in input X as the zi parameter


if __name__ == "__main__":
    # read the  raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    in_fname = "../data/test1.raw"
    iq_data = np.fromfile(in_fname, dtype="float32")
    # iq_data = (iq_data - 128.0) / 128.0
    print('Read raw RF data from "' + in_fname + '" in uint8 format')

    # coefficients for the front-end low-pass filter
    rf_coeff = our_fir(rf_Fc, rf_Fs, rf_taps)

    # coefficients for the filter to extract mono audio
    audio_coeff = our_fir(audio_Fc, 5 * audio_Fs, audio_taps)

    # set up drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace=1.0)

    # select a block_size that is in KB and
    # a multiple of decimation factors
    block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps - 1)
    state_q_lpf_100k = np.zeros(rf_taps - 1)
    state_phase = [0.0, 0.0]
    audio_zi = [0] * (audio_taps - 1)
    # add state as needed for the mono channel filter

    # audio buffer that stores all the audio blocks
    audio_data = np.array([])  # to be updated by you during in-lab

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    # while (block_count+1)*block_size < len(iq_data):
    while block_count <= 20:
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit

        print("Processing cd block " + str(block_count))
        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = our_lfilter(
            rf_coeff,
            1.0,
            iq_data[(block_count) * block_size : (block_count + 1) * block_size : 2],
            rf_taps,
            zi=state_i_lpf_100k,
        )
        q_filt, state_q_lpf_100k = our_lfilter(
            rf_coeff,
            1.0,
            iq_data[
                (block_count) * block_size + 1 : (block_count + 1) * block_size : 2
            ],
            rf_taps,
            zi=state_q_lpf_100k,
        )

        # downsample the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]
        print("----------------------------------------")
        # FM demodulator
        # fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

        fm_demod, state_phase = Demod(i_ds, q_ds, state_phase)
        # print(state_p ase)
        # extract the mono audtio data through filtering
        # audio_filt, audio_zi = our_lfilter(
        #     audio_coeff, 1.0, fm_demod, audio_taps, zi=audio_zi
        # )

        audio_filt, audio_zi = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=audio_zi)

        # downsample audio data
        audio_block = audio_filt[::audio_decim]

        # concatanete most recently processed audio_block
        # to the previous blocks stored in audio_data
        audio_data = np.concatenate((audio_data, audio_block))

        # to save runtime select the range of blocks to log iq_data
        # this includes both saving binary files as well plotting PSD
        # below we assume we want to plot for graphs for blocks 10 and 11
        if block_count == 3:
            #     print("size of x " + str(len(fm_demod)))
            #     print("zi size " + str(len(audio_zi)))
            #     print("number of taps " + str(audio_taps))

            #     print(i_filt[0:100])
            print(fm_demod[0:100])
            # PSD after FM demodulation
            # print(audio_block)
            ax0.clear()
            ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs / rf_decim) / 1e3)
            ax0.set_ylabel("PSD (dB/Hz)")
            ax0.set_xlabel("Freq (kHz)")
            ax0.set_title("Demodulated FM (block " + str(block_count) + ")")
            # output binary file name (where samples are written from Python)
            fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
            # create binary file where each sample is a 32-bit float
            fm_demod.astype("float32").tofile(fm_demod_fname)

            # PSD after extracting mono audio
            ax1.clear()
            ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs / rf_decim) / 1e3)
            ax1.set_ylabel("PSD (dB/Hz)")
            ax1.set_xlabel("Freq (kHz)")
            ax1.set_title("Extracted Mono")

            # PSD after decimating mono audio
            ax2.clear()
            ax2.psd(audio_block, NFFT=512, Fs=audio_Fs / 1e3)
            ax2.set_ylabel("PSD (dB/Hz)")
            ax2.set_xlabel("Freq (kHz)")
            ax2.set_title("Mono Audio")

            # save figure to file
            fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

        block_count += 1

    print("Finished processing the raw I/Q samples")

    # write audio data to a .wav file (assumes audio_data samples are -1 to +1)
    wavfile.write(
        "../data/fmMonoBlock.wav", int(audio_Fs), np.int16((audio_data / 2) * 32767)
    )

    # uncomment assuming you wish to show some plots
    plt.show()
