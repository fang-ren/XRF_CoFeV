# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:37:51 2016

@author: fangren
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename
from scipy.optimize import curve_fit
import imp


def func(x, *params):
    """
    create a Gaussian fitted curve according to params
    """
    y = np.zeros_like(x)
    # V alpha
    ctr1 = params[0]
    amp1 = params[1]
    wid1 = params[2]
    # V beta, refer to the paper (Smith_et_al-1974-X-Ray_Spectrometry) for the value of beta/alpha ratio
    ctr2 = params[3]
    amp2 = params[1]*0.132
    wid2 = params[4]
    # Fe alpha
    ctr3 = params[5]
    amp3 = params[6]
    wid3 = params[7]
    # Fe beta, refer to the paper (Smith_et_al-1974-X-Ray_Spectrometry) for the value of beta/alpha ratio
    ctr4 = params[8]
    amp4 = params[6]*0.136
    wid4 = params[9]
    # Co alpha
    ctr5 = params[10]
    amp5 = params[11]
    wid5 = params[12]
    # Co beta, refer to the paper (Smith_et_al-1974-X-Ray_Spectrometry) for the value of beta/alpha ratio
    ctr6 = params[13]
    amp6 = params[11]* 0.137
    wid6 = params[14]

    y = y + \
        amp1 * np.exp( -((x - ctr1)/wid1)**2) + \
        amp2 * np.exp( -((x - ctr2)/wid2)**2) + \
        amp3 * np.exp( -((x - ctr3)/wid3)**2) + \
        amp4 * np.exp( -((x - ctr4)/wid4)**2) + \
        amp5 * np.exp( -((x - ctr5)/wid5)**2) + \
        amp6 * np.exp( -((x - ctr6)/wid6)**2)
    return y



# path
path = 'mca_data\\'
base_filename = 'SampleB2_22_24x24_t30_scan1_mca.dat'

filename = path + base_filename
save_path = path + 'XRF_channels_fit\\'
if not os.path.exists(save_path):
    os.makedirs(save_path)

data = np.genfromtxt(filename, delimiter=' ')
x = np.array(range(data.shape[1]))


# # visualize a few XRF spectra for defining peak positions.
# y1 = data[1, :]
# y2 = data[100, :]
# y3 = data[300, :]
# y4 = data[440, :]
# plt.plot(x, y1)
# plt.plot(x, y2)
# plt.plot(x, y3)
# plt.plot(x, y4)


# initialize guesses and high and low bounds for fitting, refer to the excel file for peak calibration
guess = [529, 1000, 10,
         581, 10,
         684, 1000, 10,
         755, 10,
         741, 1000, 10,
         818, 10]
high = [534, 10000, 50,
         601, 100,
         704, 10000, 100,
         775, 100,
         761, 10000, 100,
         838, 100]
low = [509, 0, 0,
         561, 0,
         664, 0, 0,
         735, 0,
         721, 0, 0,
         798, 0]

# fitting starts from here
for i in range(1, data.shape[0]):
    print i
    intensity = data[i]
    try:
        popt, pcov = curve_fit(func, x, intensity, p0=guess, bounds = (low, high))
        fit = func(x, *popt)
        plt.figure(1)
        plt.plot(x, intensity)
        plt.plot(x, fit)
        ctr1 = popt[0]
        amp1 = popt[1]
        wid1 = popt[2]

        ctr2 = popt[3]
        amp2 = popt[1] * 0.132
        wid2 = popt[4]

        ctr3 = popt[5]
        amp3 = popt[6]
        wid3 = popt[7]

        ctr4 = popt[8]
        amp4 = popt[6] * 0.136
        wid4 = popt[9]

        ctr5 = popt[10]
        amp5 = popt[11]
        wid5 = popt[12]

        ctr6 = popt[13]
        amp6 = popt[11] * 0.137
        wid6 = popt[14]

        curve1 = amp1 * np.exp( -((x - ctr1)/wid1)**2)
        curve2 = amp2 * np.exp( -((x - ctr2)/wid2)**2)
        curve3 = amp3 * np.exp( -((x - ctr3)/wid3)**2)
        curve4 = amp4 * np.exp( -((x - ctr4)/wid4)**2)
        curve5 = amp5 * np.exp( -((x - ctr5)/wid5)**2)
        curve6 = amp6 * np.exp( -((x - ctr6)/wid6)**2)
        plt.plot(x, curve1)
        plt.plot(x, curve2)
        plt.plot(x, curve4)
        plt.plot(x, curve3)
        plt.plot(x, curve5)
        plt.plot(x, curve6)
        plt.xlim(380, 1000)
        print 'saving', save_path + base_filename[:-13] + str(i) + '_XRF_fit'
        plt.savefig(save_path + base_filename[:-13] + str(i) + '_XRF_fit')
        plt.close()
        result = [[amp1, amp2, amp3, amp4, amp5, amp6], [ctr1, ctr2, ctr3, ctr4, ctr5, ctr6], [wid1, wid2, wid3, wid4, wid5, wid6]]
        np.savetxt(save_path + base_filename[:-13] + str(i) + '_XRF_fit.csv', result, delimiter=",")
    except RuntimeError:
        print "Failed to fit", i+1
        print "used the previous peak information"
        np.savetxt(save_path + base_filename[:-13] + str(i) + '_XRF_fit.csv', result, delimiter=",")

