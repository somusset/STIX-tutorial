#!/usr/bin/python3
"""
 Script to plot STIX preview images  with matplotlib

 Created at 2022-05-25T08:04:31.849791 by STIX data center online image reconstruction tool
 Contact: hualin.xiao@fhnw.ch

"""
import wget
import sunpy
import sunpy.map
from astropy import units as u
from matplotlib import pyplot as plt
#%matplotlib notebook # for jupyter notebook

# define more units
u.add_enabled_units(
    [u.def_unit("arcsecs", 1 * u.arcsec),
     u.def_unit("meters", 1 * u.m)])

download_location = '.'

#download image fits files from STIX data center
fits_urls = ['https://datacenter.stix.i4ds.net/image-archive/2111010024/stix_ql_image_sci_5316_uid_2111010024_4-10keV_20211101T013654_1938_bp_map.fits', 'https://datacenter.stix.i4ds.net/image-archive/2111010024/stix_ql_image_sci_5316_uid_2111010024_4-10keV_20211101T013654_1938_vis_fwdfit_map.fits', 'https://datacenter.stix.i4ds.net/image-archive/2111010024/stix_ql_image_sci_5316_uid_2111010024_4-10keV_20211101T013654_1938_em_map.fits', 'https://datacenter.stix.i4ds.net/image-archive/2111010024/stix_ql_image_sci_5316_uid_2111010024_4-10keV_20211101T013654_1938_clean_map.fits', 'https://datacenter.stix.i4ds.net/image-archive/2111010024/stix_ql_image_sci_5316_uid_2111010024_4-10keV_20211101T013654_1938_full_disk_bp_map.fits']
filenames = [wget.download(url, out=download_location) for url in fits_urls]

maps = sunpy.map.Map(filenames)

print(f"total number of maps: {len(maps)}")

#plotting images using sunpy.map
for m in maps:
    plt.figure()
    m.plot(cmap="std_gamma_2")
    m.draw_grid(color='w', ls='--', grid_spacing=10 * u.deg)
    m.draw_limb(color='w')
plt.show()