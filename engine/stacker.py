from PIL import Image
import numpy as np
from teslaswarm.settings import STATIC_OS_PATH
import matplotlib.pyplot as plt


def stack_images(size, img1, img2, img3, img4):
    FIELD = np.zeros((size*2, size*2, 4))
    print('FIELD shape:', FIELD.shape)
    fourth1 = [
        [0, size],
        [0, size]
    ]
    fourth2 = [
        [size, size*2],
        [0, size]
    ]
    fourth3 = [
        [0, size],
        [size, size*2]
    ]
    fourth4 = [
        [size, size*2],
        [size, size*2],

    ]
    fourths = [fourth1, fourth2, fourth3, fourth4]

    for f, img in enumerate([img1, img2, img3, img4]):
        fourth = fourths[f]
        #   print(fourth[0][0],fourth[0][1], fourth[1][0],fourth[1][1])
        img = np.array(img)
        #   print(img.shape)
        FIELD[fourth[0][0]:fourth[0][1], fourth[1][0]:fourth[1][1]] = img
    return Image.fromarray((FIELD).astype(np.uint8))

def single_image(img):
    img = np.array(img)
    return Image.fromarray((img).astype(np.uint8))

#img1 = Image.open(STATIC_OS_PATH + '/media/images/2.png')
#out_image = stack_images(2000, img1, img1, img1, img1)
#out_image.save(STATIC_OS_PATH + '/media/images/3.png')