## Chromatic aberration (CA) correction software (C/C++)

#### Input parameters

Depending on the number of paramters, the program performs different tasks:

###### Three input arguments 
* `fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt`  
* Based on raw calibration image, the program estimates correction polynomials for red and bleu channels and saves them to chosen txt files  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt`  

###### Six input arguments  
* `fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm`
* Read a raw image that is needed to be corrected, reads correction polynomials and performs the correction of the image; three corrected channels are saved separately
* **EXAMPLE**: `data/_MG_7628.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`

###### More than six input arguments  
* `fname_raw_calib.pgm fname_raw_calib_red_corr.pgm fname_raw_calib_green_corr.pgm fname_raw_calib_blue_corr.pgm fname_raw_calib_keyp_dist.txt fname_raw_calib_keyp_corr.txt [fname_img_n.pgm fname_img_n_red_corr.pgm fname_img_n_green_corr.pgm fname_img_n_blue_corr.pgm, ...]`
* Runs all the circuit: first it builds the correction polynomials based on calibration raw image, it corrects the calibration image and then it corrects all the images that were taken with the same camera settings. In txt files, it saves the geometrical coordinates of the calibration image keypoints in all the channels - before and after the correction. Those txt files can be used later for visualization tests in Matlab (i.e., misalignment field and histograms). The correction polynomials are not saved, for this, use point *three-input-arguments run*.
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_R_corr.pgm data/_MG_7626_G_corr.pgm data/_MG_7626_B_corr.pgm data/_MG_7626_dist.txt data/_MG_7626_corr.txt data/_MG_7628.pgm data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`

