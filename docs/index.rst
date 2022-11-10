.. HIIdentify documentation master file, created by
   sphinx-quickstart on Thu Jun 30 09:24:35 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


..
   .. toctree::
	  :maxdepth: 2
	  :caption: Contents:



   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. image:: https://raw.githubusercontent.com/BethanEaseman/HIIdentify/master/Images/HIIdentify-logo.png
   :height: 165
   :width: 200
   :alt: HIIdentify logo


HIIdentify
==========
|

.. image:: https://readthedocs.org/projects/hiidentify/badge/?version=lateststyle=plastic
   :target: https://hiidentify.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://img.shields.io/github/last-commit/BethanEaseman/HIIdentify?style=plastic
   :target: https://img.shields.io/github/last-commit/BethanEaseman/HIIdentify?style=plastic
   :alt: Last commit


.. image:: https://img.shields.io/pypi/v/HIIdentify?style=plastic
   :target: https://img.shields.io/pypi/v/HIIdentify?style=plastic
   :alt: Version

.. image:: https://img.shields.io/badge/license-%20%20GNU%20GPLv3%20-green?style=plastic
   :target: https://img.shields.io/badge/license-%20%20GNU%20GPLv3%20-green?style=plastic
   :alt: License


Welcome to HIIdentify! This code identifies HII regions within a galaxy, using a map of the H :math:`\alpha` emssion line flux.

*Please note, HIIdentify is under active development - any contributions and / or feedback would be very welcome*.

``HIIdentify`` works by identifying the brightest pixels within the image, then growing the region to include the surrounding pixels with fluxes greater than the specified background flux, up to a maximum size. Where regions merge, the distance from the merging pixels to the peaks of the two regions are considered, and the pixel is assigned to the region with the closest peak.

In the below example map (*left*), the flux of the H :math:`\alpha` emission line can be seen, with the highest flux regions show in yellow, and lowest flux regions in purple. The regions identfied by ``HIIdentify`` can be seen as the red outlines. Here it can be seen that the regions are not restricted to being a particular shape, and that all regions with a peak flux above a given limit have been identified.

In the right-hand image, the segmentation map returned by ``HIIdentify`` can be seen. A 2D map is returned, with all pixels corresponding to a particular HII region set to the ID number of the region. This allows the segmentation map to be used to mask out regions of maps of other parameters, such as line fluxes or metallicity maps, pertaining to the selected HII region.


.. image:: https://raw.githubusercontent.com/BethanEaseman/HIIdentify/master/Images/NGC1483_ha_regionoutline_segmentationmap.png
   :alt: HIIdentify examples
   :align: center


=====================
Installing HIIdentify
=====================

``HIIdentify`` can be found on `PyPI <https://pypi.org/project/HIIdentify/>`_, and installed using `pip <https://pip.pypa.io/en/stable/>`_, by running::

    pip install HIIdentify

================
Using HIIdentify
================

The steps taken by ``HIIdentify`` can be seen below - pixels with flux values above the given max flux are first discarded from the map, then the brightest pixel is identified. Subject to this pixel meeting certain criteria, such as having fewer than 4 non-finite values in the surrounding 8, and the median of the surrounding pixels being greater than 30% of the brightness of the central pixel, then this pixel is considered to be the centre of a HII region.

The code then iterates outwards in roughly circular annuli, with pixels in these annuli being added to the region if they have a flux value > the given background flux. There is an additional criteria of requiring each considered pixel to have a flux value lower than the maximum flux of the previous annulus - this is because the flux from HII regions should decrease radially, and any pixels in subsequent annuli with higher flux values are therefore unlikely to represent gas belonging to the same HII region.

If the considered pixel has already been assigned to a different region, then the pixel is assigned to the region with the closest peak. If the pixel is isolated from the rest of the region, i.e. is is not adjacent to any pixels in the region, then it is rejected.

This outward iteration stops when >80% of the pixels within the region have been considered and rejected. If at least as many pixels as the given minimum have been added to the region, then it is considered to be a detection of a HII region, and the corresponding pixels of the segmentation map are set to the Region ID, and the Region ID is incremented by 1.

The next-brightest pixel is then considered, and the process repeated until the minimum flux limit is reached.

A .fits file is produced, and saved at <target directory>/<object name>_HII_segmentation_map.fits . Included in this file are 4 extensions, the first being a map of the Region IDs, like the above example. The second extension contains a map showing each pixel's distance to the peak of it's region. nB at this point the distance maps do not account for a galaxy's inclination.

An additional .fits file is saved at <target directory>/<object name>_HII_trackermaps.fits, with extensions of 'Peaks tracker', and 'Pixels tracker'. As the names suggest, these maps keep track of all pixels considered, either being considered as the peak of a HII region, or for being part of an identified region. The values within these maps are as follows:


.. csv-table:: Peaks tracker map
   :header: "Value", "Meaning"

	-5, "Pixel not considered"
	1, "Surrounded by >4 non-finite values - just noise."
    2, "Adjacent to an already-identified peak, or pixel with higher flux"
    3, "Median of surrounding pixels < 0.1 * peak_flux"
    4, "Separation from other region < given limit"
    5, "Pixel successfully identified as the peak of a HII region."

.. csv-table:: Pixel tracker map
   :header: "Value", "Meaning"
    -5, "Pixel not considered"
    1, "Pixel rejected - below flux threshold"
    2, "Pixel assigned to region"
	3, "Pixel has moved region"
    4, "Pixel considered, but too few pixels in region to meet threshold."
    5, "Pixel rejected - isolated / too far away from others."


Further documentation and examples will be added in due course, in the meantime, if you would like to use ``HIIdentify``, and have any questions or are interested in learning more, please email Bethan Easeman (be329@bath.ac.uk) and I'll be glad to talk it through with you!

======================
How to cite HIIdentify
======================


If you use ``HIIdentify`` as part of your work, please cite Easeman et al. (2022) *in prep*.
