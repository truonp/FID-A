%Orientation information - from FSLMRS spec2nii - **PT**2025
function [imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness, dim_swapped] = twix_orientation(twix_header,force_svs)

if nargin <2
    force_svs=0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dSag') && ~force_svs %for slice selective spectroscopy (MRSI)
    NormaldSag=twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sNormal,'dSag')
    NormaldSag=twix_header.MeasYaps.sSpecPara.sVoI.sNormal.dSag;
else
    NormaldSag=0.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dCor') && ~force_svs
    NormaldCor=twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sNormal,'dCor')
    NormaldCor=twix_header.MeasYaps.sSpecPara.sVoI.sNormal.dCor;
else
    NormaldCor=0.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dTra') && ~force_svs
    NormaldTra=twix_header.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sNormal,'dTra')
    NormaldTra=twix_header.MeasYaps.sSpecPara.sVoI.sNormal.dTra;
else
    NormaldTra=0.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'dInPlaneRot') && ~force_svs
    inplaneRotation=twix_header.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI,'dInPlaneRot')
    inplaneRotation=twix_header.MeasYaps.sSpecPara.sVoI.dInPlaneRot;
else
    inplaneRotation=0.0;
end

TwixSliceNormal = [NormaldSag, NormaldCor, NormaldTra];
% If all zeros, make a 'normal' orientation (e.g. for unlocalised data)
if ~any(TwixSliceNormal)
    TwixSliceNormal(1) = TwixSliceNormal(1) + 1.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'dReadoutFOV') && ~force_svs
    RoFoV=twix_header.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
    PeFoV=twix_header.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI,'dReadoutFOV')
    RoFoV=twix_header.MeasYaps.sSpecPara.sVoI.dReadoutFOV;
    PeFoV=twix_header.MeasYaps.sSpecPara.sVoI.dPhaseFOV;
else
    RoFoV=10000.0;
    PeFoV=10000.0;
end
if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'dThickness') && ~force_svs
    sliceThickness=twix_header.MeasYaps.sSliceArray.asSlice{1}.dThickness;
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI,'dThickness')
    sliceThickness=twix_header.MeasYaps.sSpecPara.sVoI.dThickness;
else
    sliceThickness=10000.0;
end

%Position info (including table position)
if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'sPosition') && ~force_svs
    PosdSag=0.0;
    if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dSag')
        PosdSag=twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag;
    end
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sPosition,'dSag')
    PosdSag=twix_header.MeasYaps.sSpecPara.sVoI.sPosition.dSag;
else
    PosdSag=0.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'sPosition') && ~force_svs
    PosdCor=0.0;
    if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dCor')
        PosdCor=twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor;
    end
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sPosition,'dCor')
    PosdCor=twix_header.MeasYaps.sSpecPara.sVoI.sPosition.dCor;
else
    PosdCor=0.0;
end

if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1},'sPosition') && ~force_svs
    PosdTra=0.0;
    if isfield(twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dTra')
        PosdTra=twix_header.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra;
    end
elseif isfield(twix_header.MeasYaps.sSpecPara.sVoI.sPosition,'dTra')
    PosdTra=twix_header.MeasYaps.sSpecPara.sVoI.sPosition.dTra;
else
    PosdTra=0.0;
end

if isfield(twix_header.MeasYaps,'lScanRegionPosSag')
    PosdSag = PosdSag + twix_header.MeasYaps.lScanRegionPosSag;
end
if isfield(twix_header.MeasYaps,'lScanRegionPosCor')
    PosdCor = PosdCor + twix_header.MeasYaps.lScanRegionPosCor;
end
if isfield(twix_header.MeasYaps,'lScanRegionPosTra')
    PosdTra = PosdTra + twix_header.MeasYaps.lScanRegionPosTra;
end

if ~force_svs
    %%MRSI - basically look in the twix header to see the number of voxels if
    %%more than one 1 MRSI, none = FID, 1 is svs

    data_size_pe = 1;
    if twix_header.Meas.lFinalMatrixSizePhase
        data_size_pe = twix_header.Meas.lFinalMatrixSizePhase;
    end
    data_size_ro = 1;
    if twix_header.Meas.lFinalMatrixSizeRead
        data_size_ro = twix_header.Meas.lFinalMatrixSizeRead;
    end
    data_size_sl = 1;
    if twix_header.Meas.lFinalMatrixSizeSlice
        data_size_sl = twix_header.Meas.lFinalMatrixSizeSlice;
    end
    base_pos_info=double([PosdSag,PosdCor,PosdTra]);
    [imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness, dim_swapped]=CSIOrientations(TwixSliceNormal,inplaneRotation,PeFoV,RoFoV,sliceThickness,data_size_pe,data_size_ro,data_size_sl,base_pos_info);
else
    %%not MRSI
    half_shift=0;
    [dColVec_vector,dRowVec_vector]=calc_prs(TwixSliceNormal,inplaneRotation,0);
    imageOrientationPatient = cat(1, dRowVec_vector(:).', dColVec_vector(:).');

    pixelSpacing = [PeFoV, RoFoV];  % Note: MATLAB uses row vectors by default

    imagePositionPatient = double([PosdSag, PosdCor, PosdTra]);  % ensure double precision

    dim_swapped = 0;

    xyzMM=[pixelSpacing,sliceThickness];

    dcm_to_nifti_orientation(imageOrientationPatient, imagePositionPatient, xyzMM,[1,1,1]);

    % 3) If required apply the half-voxel shift in the first two dimensions
    if half_shift
        Q44=apply_half_voxel_shift(Q44);
    end
end
end



function [imageOrientationPatient,imagePositionPatient, pixelSpacing,fov_sl,dim_swapped]=CSIOrientations(slice_normal, ip_rot,fov_pe,fov_ro,fov_sl, n_pe,n_ro,n_sl,base_pos)
SAGITTAL   = 0;
CORONAL    = 1;
TRANSVERSE = 2;

mo_case=class_ori(slice_normal(1),slice_normal(2),slice_normal(3));
[dColVec_vector,dRowVec_vector]=calc_prs(slice_normal,ip_rot);

if n_sl>1
    %3D MRSI
    base_pos = base_pos - slice_normal*(fov_sl/2-fov_sl/n_sl/2);
    fov_sl = fov_sl/n_sl;
end

if mo_case == SAGITTAL
    dRowVec_vector= dRowVec_vector*-1;
    pixelSpacing = [fov_ro/n_ro, fov_pe/n_pe];
    tmp=dRowVec_vector;
    dRowVec_vector=dColVec_vector;
    dColVec_vector=tmp;
    imagePositionPatient = base_pos - (dRowVec_vector * fov_pe / 2) - (dColVec_vector * fov_ro / 2)';
    dim_swapped=0;
elseif mo_case == CORONAL
    pixelSpacing = [fov_ro/n_ro, fov_pe/n_pe];
    tmp=dRowVec_vector;
    dRowVec_vector=dColVec_vector;
    dColVec_vector=tmp;
    imagePositionPatient = base_pos - (dRowVec_vector * fov_pe / 2) - (dColVec_vector * fov_ro / 2)';
    dim_swapped=0;
elseif mo_case == TRANSVERSE
    %Mirror along Row/readout direction: VB=LIN, VE+=SEG
    dRowVec_vector= dRowVec_vector*-1;
    pixelSpacing = [fov_pe/n_pe, fov_ro/n_ro];
    %Swap 1st and 2nd dimensions, for BG swap LIN and PHS, for VE+ swap SEG
    %and LIN dimensions
    dim_swapped=1;
    imagePositionPatient = base_pos - (dColVec_vector * fov_ro / 2)' - (dRowVec_vector * fov_pe / 2);
end
imageOrientationPatient = [dRowVec_vector; dColVec_vector'];
end
