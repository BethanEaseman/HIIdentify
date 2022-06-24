"""
Unit testing for the HIIdentify package.
"""
import unittest

import numpy as np
from astropy.io import fits

import create_test_images as cti
from ..HIIdentify.HIIdentify import identify_HII_regions


def read_in_fits(fname):
    """
    Convenenience function for reading in test maps
    """
    with fits.open(fname) as tmp:
        _map = tmp[0].data # pylint: disable=no-member
        _header = tmp[0].header # pylint: disable=no-member

    return _map, _header

class HIIdentifyTestCase(unittest.TestCase):
    """
    Test for the HIIdentify module.

    Checks that regions are correctly identified in test cases.
    """


    # read in example images (setUp())

    # For each of the example images: (?)

    # Use HIIdentify to identify regions in the image
    # Assert that regions found, and that conditions are met:
    # - If bkg_flux used, check no pixels with flux <= bkg are selected
    # - Check no peaks above flux_ulim are considered
    # - Check stops at flux_llim
    # - Check min_pixels condition is obeyed
    # - Check how many peaks identified for each image / correct no. identified
    # - Check that two regions are sensibly separated somehow if merged region used?
    # - For noisy test image, check that noisy pixels properly dealt with
    # - For noisy image with no peak, check this is dealt with correctly.



    def setUp(self):
        tdir = 'images_for_testing/'
        try:
            self.isr, self.isr_head = read_in_fits(f"{tdir}isolated_region.fits")
            self.nir, self.nir_head = read_in_fits(f"{tdir}noisy_isolated_region.fits")
            self.fir,  self.fir_head = read_in_fits(f"{tdir}faint_isolated_region.fits")
            self.mur, self.mur_head  = read_in_fits(f"{tdir}multiple_regions.fits")
            self.mer, self.mer_head  = read_in_fits(f"{tdir}merged_regions.fits")
            self.jun, self.jun_head  = read_in_fits(f"{tdir}just_noise.fits")
        except FileNotFoundError:
            cti.isolated_region()
            cti.noisy_isolated_region()
            cti.faint_isolated_region()
            cti.multiple_regions()
            cti.merged_regions()
            cti.just_noise()



    def test_isolated_region(self):
        """
        Test for an isolated region.
        """

        ideal_seg_map = np.array(self.isr[self.isr>0.5], dtype=int)

        identify_HII_regions(self.isr, self.isr_head,flux_llim=1, #pylint: disable=invalid-name
                             flux_ulim=20, z=1, obj_name='isolated_region', bkg_flux=0.5,
                             max_radius_kpc=None, max_radius_arcsec=None,
                             min_pixels=5, verbose=False, tdir='test_results')













if __name__ == "__main__":
    unittest.main()
