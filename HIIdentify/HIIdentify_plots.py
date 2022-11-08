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
import matplotlib.colors as colors
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







def draw_region_outlines(sel, x0=0, y0=0, x1 = None, y1= None):
    """
    Snippet adapted from https://stackoverflow.com/questions/24539296/outline-a-region-in-a-graph

    Takes in a boolean array of the regions, and returns arrays to draw on outlines

    sel: boolean array
        Mask of which regions to outline

    x0, y0: float
        Minimum x and y values (Default = 0)
    x1, y1: float
        Maximum x and y values (Default = size of array, works for maps)


    """

    if x1 is None:
        x1 = sel.shape[1]
    if y1 is None:
        y1= sel.shape[0]


    # a vertical line segment is needed, when the pixels next to each other horizontally
    #   belong to diffferent groups (one is part of the mask, the other isn't)
    # after this ver_seg has two arrays, one for row coordinates, the other for column coordinates
    ver_seg = np.where(sel[:,1:] != sel[:,:-1])

    # the same is repeated for horizontal segments
    hor_seg = np.where(sel[1:,:] != sel[:-1,:])

    # if we have a horizontal segment at 7,2, it means that it must be drawn between pixels
    #   (2,7) and (2,8), i.e. from (2,8)..(3,8)
    # in order to draw a discountinuous line, we add Nones in between segments
    line = []
    for pix in zip(*hor_seg):
        line.append((pix[1], pix[0]+0.5))
        line.append((pix[1]+0.5, pix[0]+0.5))
        line.append((np.nan,np.nan))

    # and the same for vertical segments
    for pix in zip(*ver_seg):
        line.append((pix[1]+0.5, pix[0]))
        line.append((pix[1]+0.5, pix[0]+0.5))
        line.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(line)

    # now we need to know something about the image which is shown
    # with this information we can rescale our points

    segments[:,0] = x0 + (x1-x0) * segments[:,0] / sel.shape[1]
    segments[:,1] = y0 + (y1-y0) * segments[:,1] / sel.shape[0]

    # Plot these segments using:
    # plt.plot(segments[:,0], segments[:,1], color='lightcoral', linewidth=0.5)

    return segments



def map_region_outline(galaxy_map, seg_map_fname, vmin=None, vmax=None, log_scaling=True,
                       origin='lower', cbar_label='', linewidth=2, tdir='./', galaxy_name='',
                       parameter_name=''):
    """
    Produces a plot of the galaxy map, with the outlines of the identified HII regions
    overlaid.

    galaxy_map: 2d array
        Map to be plotted, e.g. map of Ha flux
    seg_map_fname: str
        Path to segmentation map produced by HIIdentify
    vmin, vmax: floats
        Set to define the min and max of the colour scaling. Default = None.
    log_scaling: bool
        If True, the colour of the map will use log scaling, otherwise linear scaling will
        be used.
    origin: str
        Passed to plt.imshow. Default: 'lower' is used to account for fits.open flipping the
        map.
    cbar_label: str
        If given, label for the colourbar
    linewidth:float
        Linewidth for the region outlines. Default = 2.
    tdir: str
        Target directory for saving the file. Default='./'.
    galaxy_name: str
        Used for naming the output file. Default = ''.
    parameter_name: str
        Used for naming the output file. Default = ''.

    """

    fig,axes = plt.subplots(figsize=(17,15))
    if log_scaling:
        norm_kw = {'norm': colors.LogNorm(vmin=vmin, vmax=vmax)}
    else:
        norm_kw = {'vmin':vmin, 'vmax':vmax}
    plt.imshow(galaxy_map, origin=origin, **norm_kw)
    # origin = 'lower' accounts for the fact that fits.open flips the map
    cbar = plt.colorbar(shrink=0.8, extend='both')
    if cbar_label != '':
        cbar.set_label(r'H$\alpha$ flux', rotation=270, size=32, labelpad=25)
    axes.axis('off') #remove x and y ticks & labels for the map
    fig.tight_layout()



    #Draw outlines around regions
    with fits.open(seg_map_fname) as tmp:
        region_IDs = tmp['Region_IDs'].data
    regions = np.unique(region_IDs[region_IDs > 0])
    for r in regions:
        sel = region_IDs == r
        segments = draw_region_outlines(sel, x1 = sel.shape[1], y1= sel.shape[0])
        plt.plot(segments[:,0], segments[:,1], color='lightcoral', linewidth=linewidth)


    fname = tdir + galaxy_name + parameter_name + '_regionoutlines.png'

    fig.savefig(fname, bbox_inches="tight")

    print("\nSaved ", fname, "\n")
