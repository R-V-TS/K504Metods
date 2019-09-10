import ctypes
import numpy as np
import glob
from PIL import Image

libfile = glob.glob('build/lib.linux-x86_64-3.7/DCT_filter.cpython-37m-x86_64-linux-gnu.so')[0]
mylib = ctypes.CDLL(libfile)

mylib.DCT_filt.restype = None
mylib.DCT_filt.argtypes = [np.ctypeslib.ndpointer(dtype=np.double), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double]

def Filt_Image_via_DCT(image, window_size = 8, threshold = 0.012):
    sa = np.array(image.shape)

    if(len(sa) == 2):
        sa = np.append(sa, 1)

    print(sa)

    image_pixel = np.zeros((sa[0], sa[1] * sa[2]), dtype=ctypes.c_double)

    print(image_pixel.shape)

    for i in range(0, sa[0]):
        t = 0
        for j in range(0, (sa[1]*sa[2]), 3):
            #print(str(i) + " " + str(j) + " " + str(t))
            if(sa[2] > 1):
                for z in range(0, sa[2]):
                    image_pixel[i][j+z] = image[i][t][z]
            else:
                image_pixel[i][j] = image[i][j]
            t += 1

    mylib.DCT_filt(image_pixel, sa[0], sa[1], sa[2], window_size, threshold)

    image_result = np.zeros(sa, dtype=np.uint8)

    print(image_result.shape)

    for i in range(0, sa[0]):
        t = 0
        for j in range(0, sa[1]*sa[2], 3):
            for z in range(0, sa[2]):
                image_result[i][t][z] = image_pixel[i][j+z]
            t += 1

    return image_result


