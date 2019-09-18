import numpy as np
from PIL import Image
import K504Metods
import matplotlib.pyplot as plt

directory = "tampere17/color/"

result_ = np.zeros((300, 1))
mu = 0
sigma = 5

for z in range(0, 300):
    if (z + 1 < 10):
        im_name = 't00' + str(z + 1) + '.png'
    elif (z + 1 < 100):
        im_name = 't0' + str(z + 1) + '.png'
    else:
        im_name = 't' + str(z + 1) + '.png'

    im = Image.open(directory+im_name)

    im_array = np.array(im)
    noise = np.random.normal(mu, sigma, im_array.shape)
    im_array = im_array + noise
    imm = Image.fromarray(im_array.astype(np.uint8))
    imm.show()


    image = K504Metods.RAWImage(im_array.astype(np.short), K504Metods.TYPE_BGR)
    #image.printImage()
    coff = image.DCTCoefficients(16*3)

    image.ApplyDCT(8, 0.012)

    coff = image.DCTCoefficients(16*3)
    i = image.getImage(512*512*3)
    i = i.reshape((512,512,3))

    imm = Image.fromarray(i.astype(np.uint8))
    imm.show()
