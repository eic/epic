convert_fieldmap_to_covfie -i MARCO_v.7.6.2.2.11_1.7T_Magnetic_Field_Map_2024_05_02_rad_coords_cm_T.BMap.txt -o MARCO_v.7.6.2.2.11_1.7T_Magnetic_Field_Map_2024_05_02_rad_coords_cm_T.BMap.covfie --coord-type BrBz --r-min 0 --r-max 998 --r-step 2 --z-min -800 --z-max 798 --z-step 2

convert_fieldmap_to_covfie -i LumiDipoleMapping_2023_09_15_XYZ_coords_cm_T.txt -o LumiDipoleMapping_2023_09_15_XYZ_coords_cm_T.covfie --coord-type BxByBz --x-min -7.5 --x-max 7.5 --x-step 0.5 --y-min -34 --y-max 34 --y-step 2 --z-min -80 --z-max 80 --z-step 2
