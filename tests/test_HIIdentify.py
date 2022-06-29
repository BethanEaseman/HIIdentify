"""
Unit testing for the HIIdentify package.
"""
import unittest

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

import create_test_images as cti
from HIIdentify.HIIdentify import identify_HII_regions


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
    def setUp(self):
        self.isr, self.isr_head = cti.isolated_region()
        self.nir, self.nir_head = cti.noisy_isolated_region()
        self.fir, self.fir_head = cti.faint_isolated_region()
        self.mur, self.mur_head = cti.multiple_regions()
        self.mer, self.mer_head = cti.merged_regions()
        self.jun, self.jun_head = cti.just_noise()


    def test_isolated_region(self):
        """
        Test for an isolated region.
        """

        identify_HII_regions(self.isr.copy(), self.isr_head, flux_llim=1, #pylint: disable=invalid-name
                             flux_ulim=20, z=0.01, obj_name='isolated_region', bkg_flux=0.5,
                             max_radius_kpc=1, max_radius_arcsec=None,
                             min_pixels=5, verbose=False, tdir='test_results/')

        ideal_map = np.zeros_like(self.isr)
        ideal_map[self.isr > 0.5] = 1

        with fits.open('test_results/isolated_region_HII_segmentation_map.fits') as tmp:
            seg_map = tmp['Region_IDs'].data  # pylint: disable=no-member

        #Check that only one region identified
        self.assertTrue(np.nanmax(seg_map) == 1)

        #Check that > 90% of the pixels identified in the real map match the ideal map:
        num_matching = np.count_nonzero(ideal_map == seg_map)
        ideal_total = np.count_nonzero(self.isr > 0.5)
        seg_total = np.count_nonzero((seg_map > 0))


        self.assertTrue(num_matching > (0.9 * ideal_total))
        self.assertTrue(num_matching > (0.9 * seg_total))






    def test_noisy_isolated_region(self):
        """
        Test for a noisy isolated region
        """

        identify_HII_regions(self.nir.copy(), self.nir_head, flux_llim=9, # pylint: disable=invalid-name
                             flux_ulim=20, z=0.01, obj_name='noisy_isolated_region', bkg_flux=0.5,
                             max_radius_kpc=1, max_radius_arcsec=None,
                             min_pixels=5, min_distance=.1, verbose=False, tdir='test_results/')

        with fits.open('test_results/noisy_isolated_region_HII_segmentation_map.fits') as tmp:
            seg_map = tmp['Region_IDs'].data  # pylint: disable=no-member

        ideal_map = np.zeros_like(self.nir)
        ideal_map[self.nir > 0.5] = 1

        #Check that only one region identified
        self.assertTrue(np.nanmax(seg_map) == 1)

        num_matching = np.count_nonzero(ideal_map == seg_map)
        ideal_total = np.count_nonzero(self.nir > 0.5)
        seg_total = np.count_nonzero((seg_map > 0))

        print(f"Noisy isolated region: {ideal_total} pixels identified in the ideal map, \
{seg_total} in the seg map")

        #Check that the seg map doesn't identify more pixels as belonging to the region than the
        # ideal map
        self.assertTrue(seg_total < ideal_total)

        #Check that > 50% of the pixels identified in the real map match the ideal map:
        self.assertTrue(num_matching > (0.5 * ideal_total))



    def test_just_noise(self):
        """
        Test that no regions are identified in the just_noise map
        """

        identify_HII_regions(self.jun.copy(), self.jun_head, flux_llim=9, # pylint: disable=invalid-name
                             flux_ulim=20, z=0.01, obj_name='just_noise', bkg_flux=5,
                             max_radius_kpc=.5, max_radius_arcsec=None,
                             min_pixels=5, min_distance=.1, verbose=False, tdir='test_results/')


        with fits.open('test_results/just_noise_HII_segmentation_map.fits') as tmp:
            seg_map = tmp['Region_IDs'].data  # pylint: disable=no-member


        #Check that no regions identified
        self.assertTrue(np.nanmax(seg_map) == -50)




    def test_multiple_regions(self):
        """
        Test for multiple isolated regions
        """

        identify_HII_regions(self.mur.copy(), self.mur_head, flux_llim=9, # pylint: disable=invalid-name
                             flux_ulim=20, z=0.01, obj_name='multiple_regions', bkg_flux=0.5,
                             max_radius_kpc=1, max_radius_arcsec=None,
                             min_pixels=5, min_distance=.3, verbose=False, tdir='test_results/')


        ideal_map = np.zeros_like(self.mur)
        ideal_map[self.mur > 0.5] = 1

        with fits.open('test_results/multiple_regions_HII_segmentation_map.fits') as tmp:
            seg_map = tmp['Region_IDs'].data  # pylint: disable=no-member


        #Check that only two regions identified
        self.assertTrue(np.nanmax(seg_map) == 2)

        num_identified = np.count_nonzero((ideal_map > 0) & (seg_map > 0))
        ideal_total = np.count_nonzero(self.mur > 0.5)
        seg_total = np.count_nonzero((seg_map > 0))

        print(f"Multiple regions: {ideal_total} pixels identified in the ideal map, {seg_total} \
in the seg map")

        #Check that the seg map doesn't identify more pixels as belonging to the region than the
        # ideal map
        self.assertTrue(seg_total < ideal_total)

        #Check that > 80% of the pixels identified in the real map match the ideal map:
        self.assertTrue(num_identified > (0.8 * ideal_total))








if __name__ == "__main__":
    unittest.main()
