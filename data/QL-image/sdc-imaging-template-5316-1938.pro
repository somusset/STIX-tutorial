;Created at 2022-05-25T08:04:09.338008 by STIX data center online image reconstruction tool
;To run this script, sswidl and stix IDL software must be installed on your computer.  
; Download FITS files from STIX data center
wget("https://datacenter.stix.i4ds.net/download/fits/filename/solo_L1A_stix-sci-xray-l1-2110310071_20211031T111705-20211031T125305_035363_V01.fits",filename="solo_L1A_stix-sci-xray-l1-2110310071_20211031T111705-20211031T125305_035363_V01.fits") ; background file,
wget("https://datacenter.stix.i4ds.net/download/fits/filename/solo_L1A_stix-sci-xray-l1-2111010024_20211101T010049-20211101T023158_017507_V01.fits", filename="solo_L1A_stix-sci-xray-l1-2111010024_20211101T010049-20211101T023158_017507_V01.fits") ; signal file,
; Uncomment the following lines if you don't have stix_image_reconstruction.pro and stixmap2fits.pro on your local disk,
;wget("https://datacenter.stix.i4ds.net/pub/misc/stix_imaging/stx_image_reconstruction.pro", filename="stx_image_reconstruction.pro"),
;wget("https://datacenter.stix.i4ds.net/pub/misc/stix_imaging/stx_map2fits.pro", filename="stix_map2fits.pro"),
;.run stx_image_reconstruction.pro
;.run stix_map2fits.pro

sig_filename="solo_L1A_stix-sci-xray-l1-2111010024_20211101T010049-20211101T023158_017507_V01.fits"
bkg_filename="solo_L1A_stix-sci-xray-l1-2110310071_20211031T111705-20211031T125305_035363_V01.fits"


path_sci_file="./"+sig_filename
path_bkg_file="./"+bkg_filename
start_utc='2021-11-01T01:36:54.146'
end_utc='2021-11-01T01:37:54.146'

bp_elow=6 ; back-projection energy range lower limit
bp_ehigh=10
; energy range in  units of keV, used to make a back-project full image
; the result will be used to locate the source(s)

elow=4
ehigh=10
;energy range for EM, BP and forward-fit


; s/c emphemeris data, computed using SPICE kernel toolkits.
B0=2.111175596742068
L0=-1.7183308905557624
RSUN=1160.2317450654216
; apparent radius  of the sun in  units of arcsec
dsun=123679525985.72098
;distance between the sun and s/c in units of meters
roll_angle=-21.071092417641243
;Spacecraft roll angle in units of degrees. 
sun_center_x=-0.2036030698807715
sun_center_y=-0.5338719331034837
x_offset_arcsec=-sun_center_x
y_offset_arcsec=-sun_center_y
;Note that the off-pointing should be further corrected using the stix aspect solution


vis_fwdfit_source_type='circle'
;Change the source shape if necessary. The source shape can be also "ellipse" or "multi"  (multi-circle).


bp_fname="bp_map.fits" 
full_disk_bp_fname="full_disk_bp_map.fits" 
vis_fwdfit_fname= "vis_fwdfit_map.fits" 
em_fname= "em_map.fits"
clean_fname="clean_map.fits"
;Output filenames

stx_image_reconstruct, path_bkg_file, path_sci_file, $
	start_utc, end_utc, $
	elow, ehigh, $
	bp_elow, bp_ehigh, $
	full_disk_bp_fname,  $
	bp_fname, $
	vis_fwdfit_fname, vis_fwdfit_source_type, $
	em_fname, $
	clean_fname,  $
	L0, B0, RSUN, roll_angle, dsun, $
	x_offset_arcsec, y_offset_arcsec     


