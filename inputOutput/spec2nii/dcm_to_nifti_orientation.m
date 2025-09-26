function nii_orientation=dcm_to_nifti_orientation(imageOrientationPatient, imagePositionPatient, xyzMM, data_shape, half_shift)
if nargin<5
    half_shift = 0;
end

% in style of dcm2niix
% 1) calculate Q44
Q44=nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM);


% 2) calculate nifti quaternion parameters
% From github.com/rordenlab/dcm2niix/blob/
% 081c6300d0cf47088f0873cd586c9745498f637a/console/nii_dicom.cpp#L604
[~,Q44]=verify_slice_dir(Q44,data_shape,imagePositionPatient);
Q44(1:2,:) = Q44(1:2,:)*-1;

% 3) If required apply the half-voxel shift in the first two dimensions
if half_shift
    Q44 = apply_half_voxel_shift(Q44);
end

% 4) place in data class for nifti orientation parameters
nii_orientation = NIFTIOrient(Q44);

end

function Q44=nifti_dicom2mat(orient, patientPosition, xyzMM)
% DICOM values to 4x4 nifti transformation matrix.
% As per
% https://github.com/rordenlab/dcm2niix/blob/7ce33ca5fa3bb2dd4e5410bd97bcf515c9e462d9/console/nifti1_io_core.cpp"""" + ...
%
Q = zeros(3,3);
Q(1:2,:)=orient;

% Normalize rows
normQ = vecnorm(Q, 2, 2);           % Compute L2 norm for each row
normQ(normQ == 0.0) = 1.0;          % Replace zeros with ones
Q = Q ./ normQ;                     % Normalize each row

% row 3 is the cross product of rows 1 and 2
Q(3, :) = cross(Q(1, :), Q(2, :));

%Transpose
Q = Q.';

%if determinant Q is negative, if it is multiply third column by -1
if det(Q) < 0.0
    Q(:, 3) = Q(:, 3) * -1;
end

% next scale matrix
% dcm2niix reverses the pixel spacing?
% https://github.com/rordenlab/dcm2niix/blob/485c387c93bbca3b29b93403dfde211c4bc39af6/console/nii_dicom.cpp#L5403
tmp=xyzMM(1);
xyzMM(1)=xyzMM(2);
xyzMM(2)=tmp;

diagVox = diag(xyzMM);
Q=Q*diagVox;

Q44=zeros(4,4);
Q44(1:3,1:3)=Q;
Q44(1:3,4)=patientPosition;
Q44(4,4)=1;
end

function [iSL,R]=verify_slice_dir(R, dim, patPos)
%returns slice direction: 0 = sag, 1=coronal, 2=axial, -=flipped
%if single slice, we don't care about direction
if dim(3)<2
    iSL=[];
    return
end

iSL = 1; %find Z-slice direction; row with highest magnitude of 3rd column
if abs(R(2, 3)) >= abs(R(1, 3)) && abs(R(2, 3)) >= abs(R(3, 3))
    iSL = 2;
end
if abs(R(3, 3)) >= abs(R(1, 3)) && abs(R(3, 3)) >= abs(R(2, 3))
    iSL = 3;
end

pos = patPos(iSL);
x=[0.0,0.0,dim(3)-1,1.0];

pos1v=x*R.';
pos1=pos1v(iSL+1);

% flip = (pos > R(iSL, 4)) ~= (pos1 > R(iSL, 4));
flip = pos1 < R(iSL, 4);

if flip
    R(:,3) = R(:,3)*-1;
    iSL=iSL*-1;
end
end

function Q44=apply_half_voxel_shift(Q44)
v = [0.5, 0.5, 0] * Q44(1:3, 1:3)';
Q44(1:3, 4) = Q44(1:3, 4) + v';
end