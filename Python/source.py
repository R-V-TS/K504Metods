import numpy as np
from PIL import Image
import K504Metods
import matplotlib.pyplot as plt
from ImQuaolity import PSNR

directory = "../Images/tampere17/color/"

result_ = np.zeros((300, 1))
mu = 0
sigma = 5

for z in range(0, 1):
    if (z + 1 < 10):
        im_name = 't00' + str(z + 1) + '.png'
    elif (z + 1 < 100):
        im_name = 't0' + str(z + 1) + '.png'
    else:
        im_name = 't' + str(z + 1) + '.png'

    im = Image.open(directory+im_name)

    im_array = np.array(im)
    image = K504Metods.RAWImage(im_array.astype(np.uint8), K504Metods.TYPE_BGR)
    #image.printImage()
    #coff = image.DCTCoefficients(16*3)
    image.AddNoise(0, 5)
    i = image.getImage(512 * 512 * 3)
    i = i.reshape((512, 512, 3))
    im = Image.fromarray(i.astype(np.uint8))
    im.show()
    image.printImageCharacteristics()
    print(PSNR(im_array, i))

    image.ApplyDCT(8, 0.12)

    image.printImageCharacteristics()

    i = image.getImage(512 * 512 * 3)
    i = i.reshape((512, 512, 3))
    im = Image.fromarray(i.astype(np.uint8))
    im.show()