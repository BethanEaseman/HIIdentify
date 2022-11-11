
====================
How HIIdentify works
====================

.. note::
   The steps taken by ``HIIdentify`` can alternatively be viewed as a `flowchart <https://raw.githubusercontent.com/BethanEaseman/HIIdentify/master/Images/HIIdentify_flowchart.pdf>`_.



To identify the HII regions in a map of the H :math:`\alpha` emssion line flux, ``HIIdentify`` first discards any pixels with flux values above the given max flux, then the brightest pixel is identified. If this pixel fails to meet certain criteria, such as having fewer than 4 non-finite values in the surrounding 8, and the median of the surrounding pixels being greater than 30% of the brightness of the central pixel, then this pixel is rejected, and the next-brightest pixel identified.

Otherwise, the pixel is considered as the peak of a HII region, and the code then iterates outwards in roughly circular annuli, with pixels in these annuli being added to the region if they have a flux value > the given background flux. There is an additional criteria of requiring each considered pixel to have a flux value lower than the maximum flux of the previous annulus - this is because the flux from HII regions should decrease radially, and any pixels in subsequent annuli with higher flux values are therefore unlikely to represent gas belonging to the same HII region.

If the considered pixel has already been assigned to a different region, then the pixel is assigned to the region with the closest peak. If the pixel is isolated from the rest of the region, i.e. is is not adjacent to any pixels in the region, then it is rejected.

This outward iteration stops when >80% of the pixels within the region have been considered and rejected. If at least as many pixels as the given minimum have been added to the region, then it is considered to be a detection of a HII region, and the corresponding pixels of the segmentation map are set to the Region ID, and the Region ID is incremented by 1.

The next-brightest pixel is then considered, and the process repeated until the minimum flux limit is reached.

A .fits file is produced, and saved at <target directory>/<object name>_HII_segmentation_map.fits . Included in this file are 4 extensions, the first being a map of the Region IDs, like the above example. The second extension contains a map showing each pixel's distance to the peak of it's region. nB at this point the distance maps do not account for a galaxy's inclination.

An additional .fits file is saved at <target directory>/<object name>_HII_trackermaps.fits, with extensions of 'Peaks tracker', and 'Pixels tracker'. As the names suggest, these maps keep track of all pixels considered, either being considered as the peak of a HII region, or for being part of an identified region. Each pixel is assigned a value, with the meanings of these values described below:


.. csv-table:: Peaks tracker map
   :header: "Value", "Meaning"

   -5, "Pixel not considered"
   1, "Surrounded by >4 non-finite values - just noise."
   2, "Adjacent to an already-identified peak, or pixel with higher flux"
   3, "Median of surrounding pixels < 0.1 * peak_flux"
   4, "Separation from other region < given limit"
   5, "Pixel successfully identified as the peak of a HII region."






.. csv-table:: Pixels tracker map
   :header: "Value", "Meaning"

   -5, "Pixel not considered"
   1, "Pixel rejected - below flux threshold"
   2, "Pixel assigned to region"
   3, "Pixel has moved region"
   4, "Pixel considered, but too few pixels in region to meet threshold."
   5, "Pixel rejected - isolated / too far away from others."
