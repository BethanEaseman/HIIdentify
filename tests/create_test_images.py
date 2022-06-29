"""
Images are created for use with unit testing.

These simulated region images are used to ensure that when HIIDI is run on these images,
the expected results are produced. Simulated using 2D Gaussian profiles for each region.

Here, images for:
- an isolated region
- a noisy isolated region
- a faint isolated region
- two isolated regions
- two merged regions
- a noisy image with no region

are created, and saved in the directory images_for_testing.
"""

import numpy as np
from astropy.io import fits
from astropy.modeling.models import Gaussian2D



def create_regions(x_mean, y_mean, stddev, amplitude):
    """
    Creates the image, with the 2D Gaussian components at the specified positions.

    Single values for each of the parameters can be given, to produce a single, isolated region.
    Alternatively, a list of parameters can be given, to produce an image with mulitple regions.

    x_mean, y_mean: int
        The x and y positions of the centre of the Gaussian
    stddev: float
        Width of the Gaussian
    amplitude:
        Amplitude of the Gaussian.
    """

    x_idxs, y_idxs = np.meshgrid(np.linspace(-50,50,50), np.linspace(-50,50,50))

    if hasattr(x_mean, '__len__'):
        #If given a list of parameters, create an image with multiple Gaussians
        gaussian = Gaussian2D(x_mean=x_mean[0], y_mean=y_mean[0], x_stddev=stddev[0],\
                              y_stddev=stddev[0], amplitude=amplitude[0])

        for xmean, ymean, std, amp in zip(x_mean[1:], y_mean[1:], stddev[1:], amplitude[1:]):

            gaussian += Gaussian2D(x_mean=xmean, y_mean=ymean, x_stddev=std, y_stddev=std,\
                              amplitude=amp)

    else:
        gaussian = Gaussian2D(x_mean=x_mean, y_mean=y_mean, x_stddev=stddev, y_stddev=stddev,\
                              amplitude=amplitude)

    return gaussian(x_idxs, y_idxs)


def add_noise(array, scale):
    """
    Given an array, adds noise of the specified scale.
    """

    noise = np.random.normal(loc=0., scale=scale, size=np.array(array).shape)

    return array + noise



def isolated_region():
    """
    Created image of isolated region.

    tdir: str
        Target directory for the images to be saved into.
    """


    gaussian = create_regions(x_mean=0., y_mean=0., stddev=10., amplitude=10.)

    hdu = fits.PrimaryHDU(gaussian)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/isolated_region.fits", overwrite=True)

    return gaussian, hdu.header

def noisy_isolated_region():
    """
    Creates image of isolated region, with added noise.
    """

    gaussian = create_regions(x_mean=0., y_mean=0., stddev=10., amplitude=10.)
    gaussian = add_noise(gaussian, scale=0.5)

    hdu = fits.PrimaryHDU(gaussian)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/noisy_isolated_region.fits", overwrite=True)



def faint_isolated_region():
    """
    Creates image of a faint isolated region, with noise.
    """

    gaussian = create_regions(x_mean=0., y_mean=0., stddev=10., amplitude=1.)
    gaussian = add_noise(gaussian, scale=0.1)

    hdu = fits.PrimaryHDU(gaussian)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/faint_isolated_region.fits", overwrite=True)

    return gaussian, hdu.header


def multiple_regions():
    """
    Creates image containing several, separated objects
    """

    gaussian = create_regions(x_mean=[-30, 20], y_mean=[-30, 20], stddev=[15,7], amplitude=[10,15])
    gaussian = add_noise(gaussian, scale=0.3)

    hdu = fits.PrimaryHDU(gaussian)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/multiple_regions.fits", overwrite=True)

    return gaussian, hdu.header

def merged_regions():

    """
    Creates image with merged regions
    """
    gaussian = create_regions(x_mean=[-15, 15], y_mean=[-15, 15], stddev=[15,15], amplitude=[10,10])
    gaussian = add_noise(gaussian, scale=0.1)

    hdu = fits.PrimaryHDU(gaussian)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/merged_regions.fits", overwrite=True)

    return gaussian, hdu.header


def just_noise():
    """
    Creates image with just noise, and a bright pixel surrounded by NaN values.
    """

    image = np.ones((50,50)) + 5.

    image = add_noise(image, scale=0.1)

    # Add in bright pixel, surrounded by NaNs, to simulate noise that could be
    # interpreted by HIIdentify as a peak

    image[15:18, 15:18] = np.nan
    image[16,16] = 10.

    hdu = fits.PrimaryHDU(image)
    hdu.header['CD1_1'] = 4e-5
    hdu.header['CD1_2'] = 4e-5
    hdu.writeto("images_for_testing/just_noise.fits", overwrite=True)

    return image, hdu.header



if __name__ == '__main__':

    isolated_region()
    noisy_isolated_region()
    faint_isolated_region()
    multiple_regions()
    merged_regions()
    just_noise()
