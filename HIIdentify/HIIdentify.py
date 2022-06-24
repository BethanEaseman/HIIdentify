"""
Script contains functions to read in a 2D MUSE flux map and identify the bright regions,
taken to be the HII regions.

read_in_map:
- Reads in 2D MUSE flux map for emission line. Normally want to be using the H alpha flux map.

identify_HII_regions:
- Identifies bright regions within the given map, following the given input parameters. These
are assumed to be the HII regions when using a H alpha map.

=======================
Citation: Easeman et al. (2022), in prep.
"""

from os import path
import sys
import numpy as np


from MUSEspec.spectrum3d import Spectrum3d

from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo #pylint: disable=no-name-in-module




def read_in_map(fname, SNlim=None, verbose=False, BPT_cut=True, z=None, datadir='./',
galaxy_name=''):
    """
    Reads in a fits file from a given filename, and returns the data and header

    fname - str
        The relative filepath to the required map (normally H alpha)
    SNlim - float
        If a float is given, the flux values for pixels with a signal / noise of less
        than this value are set to 0. Default = None.
    verbose - bool
        If True, prints additional information. Default = False.
    =============
    Returns:
        Map, header
    """
    print("Reading in ", fname)
    if not path.isfile(fname): #check file exists
        raise Exception(f"File {fname} not found.")

    with fits.open(fname) as tmp:
        data, header = tmp['MAP'].data, tmp['MAP'].header # pylint: disable=no-member
        # disabling pylint error for no-member, as it's suggesting an issue with
        # .data and .header for some reason!
        if SNlim is not None:
            std = tmp['STD'].data # pylint: disable=no-member
            # Set flux values to 0. for pixels with S/N < the given limit.
            sel = (data/std) < SNlim
            data[sel] = 0.

            print_msg(f"{np.count_nonzero(sel)} pixels removed due to S/N lim", verbose)


    if BPT_cut:
        if z is None or galaxy_name=='':
            raise Exception("read_in_map: need to specify a redshift and galaxy name if want\
                             to apply BPTcut")
        s3d = Spectrum3d(filen = datadir+galaxy_name+'_gas_cube.fits.fz', datadir=datadir,
                        target=galaxy_name, z=z, fwhm=1)

        SF_sel = np.array(s3d.BPTcut(return_line_ratios=False))
        # Returns boolean mask, with True for SF pixels
        data[~SF_sel] = np.nan


    print_msg('Map read in successfully', verbose)
    return data, header, (data/std)


def determine_bkg_flux(linemap, HaEW, distmap, hasn, hasn_lim, \
                       percentile=75, dist_lim=1.5): #pylint: disable=invalid-name
    # getting pylint to ignore that HaEW doesn't conform to snake-case style
    """
    Determines a value for the background flux, to be given to identify_HII_regions.
    Finds the mean of the Ha flux values for pixels with HaEW > 6 and < 14, returns mean + std.

    linemap- 2D array
        2D array of flux values with cuts already applied for HaEW and BPT region.
    HaEW - 2D array
        2D array of HaEW values
    """

    sel = (HaEW>6) & (HaEW<14) & np.isfinite(linemap) & (distmap < dist_lim) & (hasn>hasn_lim)

    if percentile == 'median':
        return np.nanmedian(linemap[sel])

    return np.percentile(linemap[sel], percentile)


def identify_HII_regions(linemap, header, flux_llim, flux_ulim, z, #pylint: disable=invalid-name
                         obj_name, bkg_flux, max_radius_kpc=None, max_radius_arcsec=None,
                         min_pixels=2, verbose=False, tdir='./'):
    # getting pylint to ignore that HII doesn't conform to snake-case style
    """
    Using a Ha map, identifies HII regions based on the brightest regions within the map.

    linemap - 2D array
        2D array of flux values with cuts already applied for HaEW and BPT region.
    header - FITS header
        Header for the linemap file - should have CD1_1 and CD1_2 keywords.
    flux_llim - float
        Lower limit to consider when selecting out the brightest pixels in the map
    flux_ulim - float
        Upper limit of flux to consider - set to remove noise.
    flux_fraction - float
        Fraction of the peak flux at which to stop aggregating new pixels into the HII region.
    z - float
        Redshift of the galaxy
    obj_name - str
        Name of the object.
    bkg_flux - float
        Value for the background flux, used to define the edges of the HII regions.
    max_radius_kpc or max_radius_arcsec - float (in arcsec or kpc, respectively)
        Maximum radius to consider for the HII region - will stop aggregating new pixels
        beyond this. Default is None, but one must be given.
    min_pixels - int
        If given, sets the minimum number of pixels required to be within each region. Default = 2.
    verbose - bool
        If True, prints additional information. Default = False.
    tdir - str
        Path to target directory. Default='./'


    ==========
    Returns:
    ==========
       Nothing, but saves fits file named obj_name+"_HII_segmentation_map.fits" in tdir, with
       extensions of maps of the region IDs, distance from the peak, map of the identified peaks,
       and the Ha flux map for pixels assigned to HII regions.
       Maps to track which pixels were selected for HII regions, and which pixels moved region,
       are also produced, as well as maps for which pixels were considered as potential peaks &
       if / why they were rejected. These are saved as obj_name+"_HII_trackermaps.fits" in tdir.
    """

    ## Add in a check_inputs fn? i.e. check that the max flux isn't too small, flux_llim isn't
    #  too high... though do want everything to be quite generalised - maybe just pops up a
    # warning instead of raising an exception here?

    region_fluxes = linemap.copy()
    peak_tracker = np.zeros_like(linemap, dtype=float)-5
    pixel_tracker = np.zeros_like(linemap, dtype=float)-5
    # -5 is just to make the colours show up better when viewing with QFitsView or similar

    # Remove any pixels with flux above the limit:
    print(f"Removed {np.count_nonzero(linemap > flux_ulim)} pixels with flux values above the given\
 limit of {flux_ulim:.2e}")
    linemap[linemap > flux_ulim] = 0.
    num_bright_pixels = np.count_nonzero(linemap > flux_llim)
    print(f"{num_bright_pixels} pixels have flux values above the given flux_llim of \
{flux_llim:.2e}")


    maxdist, yindx, xindx, annulus_ind, AngD, pixsky = determine_distances(z, max_radius_arcsec,\
                                                                           max_radius_kpc,\
                                                                           header, linemap.shape)


    region_ID_map, peaks_map, dist_from_peak = np.zeros_like(linemap)-50, \
                                               np.zeros_like(linemap)-50, \
                                               np.zeros_like(linemap)
    #-50 is just so there's a clear difference between Region #1 and the background when
    # viewing the map e.g. in QFitsView :)
    region_ID = 1

    continue_searching = True
    loop_counter = 1

    while continue_searching:
        sys.stdout.write(f'\r Loop number {loop_counter}, \
{(loop_counter*100/num_bright_pixels):.1f} % of bright pixels considered')
        loop_counter += 1

        # Identify brightest pixel and it's index:
        peak_flux = np.nanmax(linemap)
        if peak_flux < flux_llim:
            print("\nAll bright pixels above flux_llim have been considered. Ending search for \
HII regions.")
            continue_searching = False
            break
        peak_idx = np.argwhere(linemap == peak_flux)
        peak_idx = peak_idx[0,0], peak_idx[0,1]
        print_msg(f"Bright pixel identified at {peak_idx}", verbose)


        # Check that the surrounding 8 pixels meet various conditions
        linemap, peak_tracker_val, reject_msg = check_surrounding_pixels(linemap, peak_idx,\
                                                peak_tracker_val=peak_tracker[peak_idx])
        peak_tracker[peak_idx] = peak_tracker_val
        if linemap[peak_idx] == -99: #if peak is rejected
            print_msg(reject_msg, verbose)
            continue


        # Have identified a peak, now want to check the surrounding pixels, and add them to the
        # HII region

        #Start by determining annuli to use
        annuli, distmap = determine_annuli(yindx,xindx,peak_idx,maxdist, annulus_ind, AngD, pixsky)

        # Want to go through each pixel within each annulus, checking if they meet the criteria to
        # be added to the HII region. Will then check if the number of pixels meeting the
        # criteria meets the minimum number of pixels specified in the function call.

        add_to_region, distance, moved_pixels = add_pixels_to_region(annuli, linemap, peak_idx,
        peak_flux, bkg_flux, pixel_tracker, distmap, dist_from_peak, verbose, region_ID_map)


        print_msg(f"{len(add_to_region)} pixels assigned to HII region \n", verbose)

        if len(add_to_region) > min_pixels:
            peaks_map[peak_idx] = region_ID
            peak_tracker[peak_idx] = 4.
            for idx, dist in zip(add_to_region, distance):
                region_ID_map[idx] = region_ID
                dist_from_peak[idx] = dist
                pixel_tracker[idx] = 2. #Pixel assigned to region
            for moved in moved_pixels:
                pixel_tracker[moved] = 3. #Pixel has moved region
            region_ID += 1
        else:
            for idx in add_to_region:
                pixel_tracker[idx] = 4. #Insufficient pixels in region


        linemap[peak_idx] = -99 # So that the np.nanmax picks up the next brightest pixel.



    print(f"\nIdentified {region_ID-1} HII regions. \
{(np.count_nonzero(region_ID_map[region_ID_map>0])*100)/len(linemap.flatten()):.1f}% of \
pixels assigned to a HII region \n")

    region_fluxes[region_ID_map < 0] = -50.

    ## Save maps

    hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(region_ID_map, name = 'Region_IDs'),\
                         fits.ImageHDU(dist_from_peak, name = 'Distance from Peak'), \
                         fits.ImageHDU(peaks_map, name = 'Peaks map'), \
                         fits.ImageHDU(region_fluxes, name = 'Region Flux Map')])
    hdul.writeto(tdir+obj_name+"_HII_segmentation_map.fits", overwrite=True)



    hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(peak_tracker, name = 'Peaks tracker'),\
                         fits.ImageHDU(pixel_tracker, name = 'Pixels tracker')])
    hdul.writeto(tdir+obj_name+"_HII_trackermaps.fits", overwrite=True)

    print("\nMaps saved at:  ", tdir+obj_name+"_HII_segmentation_map.fits   ", \
        tdir+obj_name+"_HII_trackermaps.fits")

    print("\n==============\n Peaks tracker map: \n 0: Pixel not considered \
    \n 1: Surrounded by >4 non-finite values - just noise. \
    \n 2: Adjacent to an already-identified peak, or pixel with higher flux \
    \n 3: Median of surrounding pixels < 0.1 * peak_flux \
    \n 4: Pixel successfully identified as the peak of a HII region. \
    \n==============\n Pixel tracker map: \n 0: Pixel not considered \
    \n 1: Pixel rejected - below flux threshold \
    \n 2: Pixel assigned to region \n 3: Pixel has moved region \
    \n 4: Pixel considered, but too few pixels in region to meet threshold. \
    \n 5: Pixel rejected - too far away from others. \
    \n==============\n\n")




def determine_distances(z, max_radius_arcsec, max_radius_kpc, header, linemap_shape):
    """
    Takes the redshift of a galaxy, as well as the maximum radius required for the HII region,
    and returns the maximum radial size in pixels, as well as maps of the x and y index for
    each pixel in the map.

    z: float
        Redshift of the galaxy
    max_radius_kpc or max_radius_arcsec - float (in arcsec or kpc, respectively)
        Maximum radius to consider for the HII region - will stop aggregating new pixels
        beyond this.
    header: header of FITS file
        Header of the Ha map FITS file - need to contain 'CD1_1' and 'CD1_2', used to convert
        between pixels and arcseconds.
    linemap_shape: tuple
        Shape of the Ha map
    ==========
    Returns:
    maxdist: float, maximum distance from the centre of the region, in kpc
    yindx, xindx: arrays of int - used to determine the distance from the centre of the region
    ==========
    """
    #Determine the maximum distance from the peak to consider pixels for HII region
    AngD=(cosmo.kpc_proper_per_arcmin(z).value[0]) /60
    pixsky = (header['CD1_1']**2 + header['CD1_2']**2) ** 0.5 * 3600

    # arcsec to kpc - multiply by AngD
    # pixels to arcsec - mulitply by pixsky


    if (max_radius_arcsec is None) and (max_radius_kpc is None):
        raise Exception("One of max_radius_kpc or max_radius_arcsec must be given!")
    if max_radius_kpc is not None:
        # If max radius given in kpc, convert to arcsec and pixels
        max_radius_arcsec = max_radius_kpc / AngD
        max_radius_pix = max_radius_kpc / (AngD * pixsky)
    elif max_radius_arcsec is not None:
        max_radius_kpc = max_radius_arcsec * AngD
        max_radius_pix = max_radius_arcsec / pixsky

    print(f"Maximum radius for the HII regions set to {max_radius_pix:.2f} pixels, \
{max_radius_arcsec:.3f} arcsec, {max_radius_kpc:.3f} kpc")
    yindx, xindx = np.indices(linemap_shape)[0], np.indices(linemap_shape)[1]
    #Used to produce map of distance from peak


    annulus_ind = np.cumsum(np.sort(np.append(np.arange(2, 50), np.arange(2, 50)))) #Cumulative sum
    # of array [2,2,3,3,4,4...]. Used when selecting annuli of pixels moving outward from the peak,
    # to determine whether pixels belong to HII region - gives circular-ish shape. Used because the
    # first annulus out from the centre contains 2 different distances, as does the subsequent one,
    # then 3 distances...
    annulus_ind = np.insert(annulus_ind, 0, 0) #Add in the central pixel


    return max_radius_kpc, yindx, xindx, annulus_ind, AngD, pixsky



def check_surrounding_pixels(linemap, peak_idx, peak_tracker_val):
    """
    Checks that the 8 pixels surrounding the selected 'peak' meet various conditions.

    linemap: 2D array
        2D array of flux values with cuts already applied for HaEW and BPT region.
    peak_idk: tuple
        Index of the position of the peak within the map.
    peak_tracker_val: float
        Current value of the pixel in the pixel tracker map.

    ==========
    Returns:
    ==========
    peak_tracker_val: int, used in the tracker map to show if/why peaks were rejected
    reject_peak: bool, whether to reject peak
    reject_msg: str, message describing reason for rejection
    """

    surrounding_ha = linemap[peak_idx[0]-1:peak_idx[0]+2, peak_idx[1]-1:peak_idx[1]+2].flatten()
    surrounding_ha = np.delete(surrounding_ha, 4) # remove central

    if np.count_nonzero(~np.isfinite(surrounding_ha)) > 4:
        reject_msg = "Bright pixel is surrounded by >4 non-finite values"
        peak_tracker_val = 1.
        linemap[peak_idx] = -99

    elif np.any(surrounding_ha == -99):
        reject_msg = "Adjacent to an already-identified peak, or pixel with higher flux - \
disregarding this one"
        peak_tracker_val = 2.
        linemap[peak_idx] = -99

    elif np.nanmedian(surrounding_ha) < (0.3 * linemap[peak_idx]):
        reject_msg = "Disregarding peak - median of surrounding pixels < 0.3 * peak_flux"
        peak_tracker_val = 3.
        linemap[peak_idx] = -99

    else:
        reject_msg=''
        peak_tracker_val = -5


    return linemap, peak_tracker_val, reject_msg



def determine_annuli(yindx, xindx, peak_idx, maxdist, annulus_ind, AngD, pixsky):
    """
    Determine the indices of pixels falling within each annulus, moving outwards
    from the index of the peak position.

    yindx, xindx : array of int
        Indices of each element on the linemap
    peak_idx: tuple
        Indices of the identified peak
    maxdist: float
        Maximum radial distance from the centre to consider pixels for the region.
    annulus_ind: array of int
        How many distances to consider when determining indices for pixels within
        each annulus

    ==========
    Returns:
    ==========
    annuli: nested list of indices of the pixels falling within each annulus.
    """

    distmap = ((yindx - peak_idx[0])**2 + (xindx - peak_idx[1])**2) ** 0.5
    #The radial distance of each point from the peak (in pixels)

    #Convert from pixels to kpc
    # arcsec to kpc - multiply by AngD
    # pixels to arcsec - mulitply by pixsky

    distmap = distmap * pixsky * AngD

    dists = np.sort(np.unique((distmap)))
    dists = dists[dists < maxdist ] # don't consider beyond the maximum distance
    dists = dists[1:] #ignore the central pixel

    # Determine the indices of the pixels falling in each radial annulus moving outward
    # from the peak.
    slices = annulus_ind[annulus_ind < len(dists)]
    annuli = []
    for lower,upper in zip(slices[:-1], slices[1:]):
        annulus = []
        for dist in dists[lower:upper]:
            annulus.extend(list(np.argwhere(distmap == dist)))
        annuli.append(annulus)


    return annuli, distmap



# def add_pixels_to_region(annuli, linemap, peak_idx, peak_flux, flux_fraction, pixel_tracker,
# distmap, dist_from_peak, verbose, region_ID_map):
#    """
#     Determines which pixels meet the conditions to be added to the region. moving radially
#     outwards from the peak, in annuli.

#     ==========
#     Returns:
#     ==========
#     add_to_region: list of indices of pixels to be added to the region
#     distance: distance of each pixel from the peak
#     moved_pixels: list of indices of pixels which have moved region, to be marked in the pixel
#      tracker map
#     """

#     add_to_region, distance, moved_pixels = [peak_idx],[0.], []
#     med_prev_annulus = 1e100 # Want medians to get lower as move outwards

#     for annulus in annuli:
#         num_in_annulus, num_below_lim = len(annulus), 0

#         annulus_fluxes = np.array([linemap[y,x] for y,x in annulus])
#         already_assigned = np.array([region_ID_map[y,x]>0 for y,x in annulus])

#         if np.nanmedian(annulus_fluxes[~already_assigned]) > med_prev_annulus:
#             break # Stop if median has gone down and then up again
#             #Only consider unassigned pixels to avoid issues with merged regions.


#         for y,x in annulus: # Go through each pixel in the annulus

#             # If the flux is below peak flux * flux fraction, then reject
#             if not (linemap[y,x] >= (peak_flux * flux_fraction)): #so NaNs also excluded
#                 #y,x as MUSE cubes in the form y,x,lambda when read in using fits.open().
#                 num_below_lim += 1
#                 if pixel_tracker[y,x] <= 0: #if hasn't already been assigned to region
#                     pixel_tracker[y,x] = 1.
#                 continue

#             # Otherwise...
#             if region_ID_map[y,x] == -50:
#                 # If pixel hasn't been assigned to a region already.
#                 add_to_region.append((y,x))
#                 distance.append(distmap[y,x])
#             else:
#                 print(f"Pixel ({x},{y}) has already been assigned to a region", verbose)
#                 if dist_from_peak[y,x] <= distmap[y,x]:
#                     #If closer to previous peak, leave there
#                     print_msg("Closer to previous peak, leaving where it is.", verbose)
#                 else:
#                     print_msg("Moving to new region", verbose)
#                     moved_pixels.append((y,x))
#                     add_to_region.append((y,x))
#                     distance.append(distmap[y,x])

#         if num_below_lim == num_in_annulus: # If all pixels in that annulus
#             # were below the threshold, then stop there.
#             print_msg("Finished adding pixels to region", verbose)
#             break
#         med_prev_annulus = np.nanmedian(annulus_fluxes)

#     return add_to_region, distance, moved_pixels


def add_pixels_to_region(annuli, linemap, peak_idx, peak_flux, bkg_flux, pixel_tracker, distmap,
                         dist_from_peak, verbose, region_ID_map, isolated_pixel_limit=50):
    """
    Determines which pixels meet the conditions to be added to the region. moving radially outwards
    from the peak, in annuli.

    isolated_pixel_limit: float or int
        Used to ignore any isolated pixels - if a pixel is at least the given distance away from the
        rest of the pixels in the region, then it is ignored. Default = 50 (parsecs).

    ==========
    Returns:
    ==========
    add_to_region: list of indices of pixels to be added to the region
    distance: distance of each pixel from the peak
    moved_pixels: list of indices of pixels which have moved region, to be marked in the pixel
     tracker map
    """

    add_to_region, distance, moved_pixels = [peak_idx],[0.], []
    max_flux = peak_flux # Starting value.  Want to ensure fluxes get lower as move outwards.

    for annulus in annuli:
        num_outside_lim = 0
        annulus_fluxes = []

        for y,x in annulus: # Go through each pixel in the annulus

            if not (linemap[y,x] <= max_flux) & (linemap[y,x] > bkg_flux):
                #reject spaxel
                #y,x as MUSE cubes in the form y,x,lambda when read in using fits.open().
                num_outside_lim += 1
                if pixel_tracker[y,x] <= 0: #if hasn't already been assigned to region
                    pixel_tracker[y,x] = 1.
                continue

            if distmap[y,x] > np.max(distance)+isolated_pixel_limit: #ignore isolated pixel
                num_outside_lim += 1
                if pixel_tracker[y,x] <= 0: #if hasn't already been assigned to region
                    pixel_tracker[y,x] = 5. #isolated pixel
                continue

            # Otherwise, add to the region.
            if region_ID_map[y,x] == -50:
                # If pixel hasn't been assigned to a region already.
                add_to_region.append((y,x))
                distance.append(distmap[y,x])
            else:
                print_msg(f"Pixel ({x},{y}) has already been assigned to a region", verbose)
                if dist_from_peak[y,x] <= distmap[y,x]:
                    #If closer to previous peak, leave there
                    print_msg("Closer to previous peak, leaving where it is.", verbose)
                else:
                    print_msg("Moving to new region", verbose)
                    moved_pixels.append((y,x))
                    add_to_region.append((y,x))
                    distance.append(distmap[y,x])
            annulus_fluxes.append(linemap[y,x])

        if num_outside_lim > 0.8*len(annulus): # If above 80% of pixels in that annulus
            # were outside the threshold limits, then stop there.
            print_msg("Finished adding pixels to region", verbose)
            break

        if np.nanmax(annulus_fluxes) < max_flux: #+np.std(annulus_fluxes)
            max_flux = np.nanmax(annulus_fluxes)#+np.std(annulus_fluxes)


    # Remove any isolated pixels
    to_be_added = set(add_to_region)
    for y,x in to_be_added:
        surrounding = [(y-1, x-1), (y-1, x), (y-1, x+1), (y, x-1), (y, x+1), (y+1, x-1), (y+1, x),\
                       (y+1, x+1)]
        in_to_be_added = [s for s in surrounding if s in to_be_added]
        if len(in_to_be_added) == 0: # If isolated
            add_to_region.remove((y,x))
            if pixel_tracker[y,x] <= 0: #if hasn't already been assigned to region
                pixel_tracker[y,x] = 5.
            continue

        ## check that moved pixels are closer to the flux of the new region
        # if (y,x) in moved_pixels:
        #     prev_region = linemap[y-1:y+2, x-1:x+2][region_ID_map[y-1:y+2, x-1:x+2] \
        #                                              == region_ID_map[y,x]]
        #     curr_region = [linemap[i,j] for i,j in in_to_be_added]

        #     diff_prev = np.abs(linemap[y,x] - np.nanmean(prev_region)) < np.std(prev_region)
        #     diff_curr = np.abs(linemap[y,x] - np.nanmean(curr_region)) < np.std(curr_region)

        #     # If within a standard dev of both means, leave it so it's moved to the closest peak

        #     # If the difference between the flux and the mean of the previous region is smaller
        #     # than that for the current region, then don't move it.
            # if diff_curr and (np.abs(linemap[y,x] - np.nanmean(prev_region)) < \
            #                   np.abs(linemap[y,x]- np.nanmean(curr_region))):
        #         moved_pixels.remove((y,x))
        #         add_to_region.remove((y,x))
        #         distance.remove(distmap[y,x])


    return add_to_region, distance, moved_pixels



def print_msg(msg, verbose):
    """
    If verbose is True, prints the message to the screen - just def'd to tidy up the code.
    """

    if verbose:
        print(msg)
