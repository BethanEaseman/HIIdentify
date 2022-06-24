"""
Images are created for use with unit testing.

These simulated region images are used to ensure that when HIIDI is run on these images,
the expected results are produced. Simulated using 2D Gaussian profiles for each region.

Here, images for:
- an isolated region
- a noisy isolated region
- a faint isolated region
- two merged regions
- a noisy image with no region

are created, and saved in the directory images_for_testing.
"""

import matplotlib.pyplot as plt
import numpy as np

def isolated_region(tdir):
    """
    Creates a simulated image of an isolated region.

    tdir: str
        Target directory for the images to be saved into.
    """

    fig = plt.figure()


    x_idxs, y_idxs = np.meshgrid(np.linspace(-50,50,50), np.linspace(-50,50,50))
    dist = np.sqrt(x_idxs**2+y_idxs**2)
    sigma, mean = 10, 0.

    gaussian = np.exp(-( (dist-mean)**2 / ( 2.0 * sigma**2 ) ) )

    plt.imshow(gaussian)
    plt.colorbar()

    fig.savefig(f'{tdir}test.png')
    print(f" Saved figure at {tdir}test.png")






    # noisy_isolated_region(tdir)
    # faint_isolated_region(tdir)
    # merged_regions(tdir)
    # noisy_image(tdir)


if __name__ == '__main__':
    TDIR = 'images_for_testing/' # target directory

    isolated_region(TDIR)
    # noisy_isolated_region(TDIR)
    # faint_isolated_region(TDIR)
    # merged_regions(TDIR)
    # noisy_image(TDIR)
