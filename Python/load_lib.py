import ctypes
import numpy as np
import glob
from PIL import Image

libfile = glob.glob('build/lib.linux-x86_64-3.7/DCT_frequency.cpython-37m-x86_64-linux-gnu.so')[0]
mylib = ctypes.CDLL(libfile)

mylib.DCT_frequency.restype = None
mylib.DCT_frequency.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_short), ctypes.c_int, ctypes.c_int, ctypes.c_int]

image = Image.open("../Images/baboon.png")

image_pixel = np.array(image, dtype=np.int16)
print(image_pixel)
sa = np.array(image_pixel.shape)
mylib.DCT_frequency(image_pixel, sa[1], sa[0], sa[2])