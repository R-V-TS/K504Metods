from PIL import Image
import load_lib
import numpy as np

image = Image.open("baboon.png")

image_pixel = np.array(image)
print(image_pixel.shape)

im = load_lib.Filt_Image_via_DCT(image_pixel, 2, 0.5) # max window size = 32

print(im)
image_ = Image.fromarray(im)
image_.show()