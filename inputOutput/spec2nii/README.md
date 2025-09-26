Using spec2nii in MATLAB
By: Peter Truong, SRI, 2025

I managed to port over the essential functions from spec2nii (original, python code) into MATLAB. 
I focused mainly on just getting the orientation in nifti format.
The two main functions are:
	1)nii_orientation=svs_orientation(twix_header);
	2)nii_orientation=mrsi_orientation(twix_header,num_t_pts);
	
Put these lines in the load function. I put both in io_loadspec_twix and io_CSIload_twix at the end.
	
This will output the affine matrix and the bcd,xyz coordinates that are needed for writing out nifti files with correct position/orientation.
More functions on the way, just wanted to share the base so we can get started on writing out nifti files in FID-A!

Down the line, I probably want to slim down the function a bit more, to only input the necessary header variables.
Write functions to be developed!