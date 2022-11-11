.. HIIdentify documentation master file, created by
   sphinx-quickstart on Thu Jun 30 09:24:35 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


..
   .. toctree::
	  :maxdepth: 2
	  :caption: Contents:

	  ./examples.rst
	  ./how_HIIdentify_works.rst


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


===============================
Getting started with HIIdentify
===============================

:doc:`how_HIIdentify_works`
	 A more in-depth description of the steps taken by ``HIIdentify`` to identify the HII regions within the galaxy


:doc:`examples`
	 Examples of how to use ``HIIdentify``






Further documentation and examples will be added in due course, in the meantime, if you would like to use ``HIIdentify``, and have any questions or are interested in learning more, please email Bethan Easeman (be329@bath.ac.uk) and I'll be glad to talk it through with you!

======================
How to cite HIIdentify
======================


If you use ``HIIdentify`` as part of your work, please cite Easeman et al. (2022) *in prep*.
