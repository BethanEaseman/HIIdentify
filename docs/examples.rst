Here's a quick guide on how to use HIIdentify to identify HII regions within your galaxy image.

The functons within ``HIIdentify`` are:
1. ``determine_bkg_flux`` - this function takes in the H :math:`\alpha` emission line map, the H :math:`\alpha` equivalent width map, and a map of the H :math:`\alpha` signal / noise. Selections are then made to cut out pixels which do not meet the given criteria, and the given percentile of the remaining flux data is returned to be used as the background flux level.


2. ``identify_HII_regions`` - this is the main function, which identifies the regions within your map, and produces the segmentation map.


==============
Example of use
==============

For this example, we'll use the H :math:`\alpha` map for NGC1483, as observed as part of the MUSE MAD sample. The H :math:`\alpha` flux, error, and equivalent width maps, produced during the analysis presented in Easeman et al. (2022) (*in prep*), can be accessed ***.

To start with, we need to determine our input values for ``identify_HII_regions``. Note that as the physical parameters of HII regions are not well constrained, a certain amount of trial and error may be required to produce segmentation maps which produce an acceptable segmentation map.

A tool for choosing a value for the background flux exists within ``HIIdentify``,


Once the background flux has been chosen,
