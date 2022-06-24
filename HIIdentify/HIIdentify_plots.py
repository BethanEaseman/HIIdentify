"""
Contains functions to produce plots from the HIIdentify maps.

distr_of_region_sizes - Using the HIIdentify segmentation map, plots a histogram of the
    sizes of the identified regions.

HII_Region_dist_from_centre - Produces a histogram, showing the radial distribution of identified
    HII regions, with the position of selected regions marked on with vertical lines.

=======================
Citation: Easeman et al. (2022), in prep.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


# Note: if using pylint, a number of disable=no-member commands have been added, as it errorneously
# flags fits files as having no 'data' attribute.
# If you're not using pylint, just ignore them :)

def distr_of_region_sizes(seg_map_fname, tdir='./'):
    """
    Using the HIIdentify segmentation map, plots a histogram of the sizes of the identified regions.

    seg_map_fname: str
       Path to the segmenation map.
    tdir: str
       Relative path to the target directroy for the files to be saved into.
    """


    with fits.open(seg_map_fname) as tmp:
        distmap = tmp['DISTANCE FROM PEAK'].data # pylint: disable=no-member
        regions = tmp['REGION_IDS'].data # pylint: disable=no-member


    region_sizes = []
    for region in np.unique(regions):
        region_sizes.append(np.nanmax(distmap[regions==region]))


    fig = plt.figure()
    plt.hist(region_sizes, bins=20)
    plt.xlabel('Max radius of identified regions (kpc)')
    plt.ylabel('Number of regions')
    fig.savefig(f"{tdir}size_of_identified_regions.png")

    print(f"Figure saved at: {tdir}size_of_identified_regions.png ")





def radial_distribution(datadir, galaxy_name, distmap_fname, \
                                distmap_key = 'DISTMAP (R/RE)', tdir='./',\
                                selected_regions=None):
    """
    Produces a histogram, showing the radial distribution of identified HII regions, with the
    position of selected regions marked on with vertical lines.
    """


    with fits.open(datadir+galaxy_name+"_HII_segmentation_map.fits") as tmp:
        peak_pos = tmp['Peaks map'].data # pylint: disable=no-member
    with fits.open(distmap_fname) as tmp:
        distmap = tmp[distmap_key].data # pylint: disable=no-member

    regions = np.unique(peak_pos[peak_pos>0])

    dists = np.array([distmap[peak_pos == region][0] for region in regions])

    fig = plt.figure()
    plt.hist(dists, bins=30, linewidth=0.3, edgecolor='k', color='lightskyblue')

    if selected_regions is not None:

        prop_cycle = plt.rcParams['axes.prop_cycle']
        marker_colours = prop_cycle.by_key()['color']
        while len(marker_colours) < len(selected_regions):
            marker_colours.extend(marker_colours)

        for sel_reg, colour in zip(selected_regions, marker_colours):
            sel = regions==sel_reg
            plt.axvline(dists[sel], color=colour, linewidth=2, label=str(sel_reg))


    plt.xlabel("Distance from galaxy centre (R/Re)")
    plt.ylabel("Number of HII regions")
    plt.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", title="Region ID")
    fig.savefig(tdir+"Regions_dist_from_galaxy_centre.png")
    print("\n Saved figure: %s \n"%(tdir+"Regions_dist_from_galaxy_centre.png"))
