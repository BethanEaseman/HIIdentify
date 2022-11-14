
===========================
Example of using HIIdentify
===========================


``HIIdentify`` provides a number of functions to identify the HII regions with your map, as well as functions to help determine the input values.


The functons within ``HIIdentify`` are:

``determine_bkg_flux`` - this function takes in the H :math:`\alpha` emission line map, the H :math:`\alpha` equivalent width map, and a map of the H :math:`\alpha` signal / noise. Selections are then made to cut out pixels which do not meet the given criteria, and the given percentile of the remaining flux data is returned to be used as the background flux level.

``identify_HII_regions`` - this is the main function, which identifies the regions within your map, and produces the segmentation map.


For this example, I'll use the H :math:`\alpha` map for NGC1483, as observed as part of the MUSE MAD sample.

.. code-block:: python

   from HIIdentify import HIIdentify

   # Read in your H alpha map, uncertainties on the flux values, and H alpha equivalent
   # width map and uncertainties, e.g.:
   ha_maps = fits.open('<filepath>')
   ha_flux = ha_maps['FLUX'].data
   ha_err = ha_maps['ERR'].data
   HaEW = ha_maps['HAEW'].data
   HaEW_err = ha_maps['HAEW_ERR'].data


To start with, we need to determine our input values for ``identify_HII_regions``. Note that as the physical parameters of HII regions are not well constrained, a certain amount of trial and error may be required to produce segmentation maps which produce an acceptable segmentation map.

A tool for choosing a value for the background flux exists within ``HIIdentify``, if we wish to use this, we need to give it our H :math:`\alpha` maps, as well as our required criteria for various parameters. For example, we can give the lower and upper limits on the H :math:`\alpha` equivalent widths, which default to values of 6 and 14, respectively, to remove pixels associated with star formation.

The ``determine_bkg_flux`` function returns the given percentile of the flux values which meet the various HaEW and S/N criteria as a measure of the background flux level. Here we'll get the function to return the 75th percentile:

.. code-block:: python

   bkg_flux = HIIdentify.determine_bkg_flux(linemap = ha_flux, HaEW = HaEW,\
				             HaEW_sn = HaEW/HaEW_err, HaEWsnlim=15,\
				             distmap=None, hasn = ha_flux/ha_err,\
				             hasn_lim=15, percentile = 75)


Once the background flux has been chosen, it can be given to ``identify_HII_regions``. The given H :math:`\alpha` map should first be cleaned, for example setting pixels with H :math:`\alpha` equivalent width < 14 to NaN to select out pixels associated with star formation.

The header is also required - this should have 'CD1_1' and 'CD1_2' keywords, used to convert
between pixels and arcseconds.

The upper and lower flux limits within which to consider pixels as potential peaks of HII regions must also be given. The maximum size of the regions can be specified, here we'll set it to a large value of 1 kpc, so that this is not the limiting factor and instead the region grows while the flux of the pixels is larger than the background flux. The minimum separation is given in kpc.

.. code-block:: python

   HIIdentify.identify_HII_regions(Hamap, ha_header, flux_llim=3e3, flux_ulim=1e7, z=0.004,\
				   obj_name='NGC1483', bkg_flux=bkg_flux, max_radius_kpc=1,\
				   min_pixels=10, verbose=False, tdir='./', min_separation=0.05)


In our given target directory (tdir), files will be created with the segmentation maps (<tdir><object name>_HII_segmentation_map.fits), and tracker maps for the pixels and peaks (<tdir><object name>_HII_trackermaps.fits). An explanation of the values within the tracker maps can be found in :doc:`how_HIIdentify_works`.


The H :math:`\alpha` map can be plotted with the region outline overlaid, and plots showing the distribution of region sizes - more details can be found in :doc:`producing_plots`.
