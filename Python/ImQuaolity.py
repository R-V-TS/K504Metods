import numpy as np
from math import pow, log10

def MSE(P, Q):
    s1 = P.shape
    s2 = Q.shape
    if(s1 != s2):
        print("Images shapes not equal")
        return
    MSE = np.zeros(3)
    for channel in range(0, s1[2]):
        MSE_local = 0
        for i in range(0, s1[0]):
            for j in range(0, s1[1]):
                MSE_local = MSE_local + pow(P[i][j][channel] - Q[i][j][channel], 2)
        MSE_local = MSE_local / (s1[1]*s1[0])
        MSE[channel] = MSE_local

    return MSE

def PSNR(P, Q):
    s1 = P.shape
    s2 = Q.shape
    if (s1 != s2):
        print("Images shapes not equal")
        return
    MSE_im = MSE(P, Q)
    PSNR = np.zeros(3)
    for channel in range(0, s1[2]):
        PSNR_local = 0
        max = 0
        for i in range(0, s1[0]):
            for j in range(0, s1[1]):
                if(P[i][j][channel] > max): max = P[i][j][channel]

        PSNR_local = 10 * log10(max/MSE_im[channel])
        PSNR[channel] = PSNR_local

    return PSNR
