function nii_orientation=svs_orientation(twix_header)
force_svs=1;
%Perform Orientation Calculations
%1) Calculate dicom like imageOrientationPatient, imagePositionPatient, pixelSpacing and slicethickness
[imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness] = twix_orientation(twix_header,force_svs);
%2) In the style of dcm2niix calculate the affine matrix
nii_orientation=dcm_to_nifti_orientation(imageOrientationPatient,imagePositionPatient,[pixelSpacing,sliceThickness],[1,1,1],0);
end