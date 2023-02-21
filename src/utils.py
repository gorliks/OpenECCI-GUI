import datetime
import time
from PIL import Image
import cv2
import numpy as np
import cv2
import matplotlib.pyplot as plt


def map_scrollbar_to_value(x_start : float = -5,
                           x_end   : float = +5,
                           number_of_steps : int =1000,
                           value_to_map : int = 0) -> float:
    a0 = x_start
    slope = (x_end - x_start) / number_of_steps
    return (a0 + slope * value_to_map)


def map_value_to_scrollbar(x_start : float = -5,
                           x_end   : float = +5,
                           number_of_steps : int =1000,
                           value_to_map : int = 0) -> int:
    slope = number_of_steps / (x_end - x_start)
    return int( (value_to_map-x_start) * slope )


def load_image(file_path):
    #image = Image.open(file_path).convert('L') # load image .tiff .png .jpg, convert to grayscale
    #image = np.array(image, dtype=np.float64)
    image = plt.imread(file_path)
    return image


def normalise(image):
    norm = np.linalg.norm(image)
    image_norm = image / norm  # normalized matrix
    return image_norm

def calculate_difference(image1, image2):
    return normalise(image1) - normalise(image2)


def current_timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime("%y%m%d.%H%M%S")


def enhance_contrast(image, clipLimit=1.0, tileGridSize=8):
    # if image.dtype == 'float64' or 'uint64':
    #     pass

    # the images needs to be 8bit, otherwise the algorithm does not work TODO fix 8-bit to any-bit
    try:
        image = image / image.max()
    except:
        image = image.data / image.data.max()
    image = (image * 2 ** 8).astype('uint8')

    tileGridSize = int(tileGridSize)

    clahe = cv2.createCLAHE(clipLimit=clipLimit,
                            tileGridSize=(tileGridSize, tileGridSize))
    image = clahe.apply(image)
    return image


def equalise_histogram(image, bitdepth=8):
    try:
        _image = image.data
    except:
        _image = image

    # _image = _image / _image.max()
    # _image = (_image * 2 ** bitdepth).astype("uint8")

    _image = cv2.equalizeHist(_image)
    return _image


def resize(image, size=(200,200)):
    return cv2.resize(image, size)