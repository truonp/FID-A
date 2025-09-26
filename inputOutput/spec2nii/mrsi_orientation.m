function nii_orientation=mrsi_orientation(twix_header,sz_t)
force_svs=0;

data_size = [];
if twix_header.Meas.lFinalMatrixSizePhase
    data_size(end+1) = int32(twix_header.Meas.lFinalMatrixSizePhase);
else
    data_size(end+1) = int32(1);
end
if twix_header.Meas.lFinalMatrixSizeRead
    data_size(end+1) = int32(twix_header.Meas.lFinalMatrixSizeRead);
else
    data_size(end+1) = int32(1);
end
if twix_header.Meas.lFinalMatrixSizeSlice
    data_size(end+1) = int32(twix_header.Meas.lFinalMatrixSizeSlice);
else
    data_size(end+1) = int32(1);
end
data_size(end+1)=sz_t;

%Perform Orientation Calculations
%1) Calculate dicom like imageOrientationPatient, imagePositionPatient, pixelSpacing and slicethickness
[imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness, dim_swapped] = twix_orientation(twix_header,force_svs);
%2) In the style of dcm2niix calculate the affine matrix
nii_orientation=dcm_to_nifti_orientation(imageOrientationPatient,imagePositionPatient,[pixelSpacing,sliceThickness],data_size(1:3),1);
end