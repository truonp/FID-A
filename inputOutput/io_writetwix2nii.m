function out=io_writetwix2nii(in)
%Parsing FID-A format to NIfTI-MRS format as described in:
%`Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher
%BJ, Wilson M. NIfTI-MRS: A standard data format for magnetic resonance
%spectroscopy. Magn Reson Med. 2022. doi: 10.1002/mrm.29418.` & https://www.nitrc.org/forum/message.php?msg_id=3738


% C Header Struct
% /*! \struct nifti_2_header
%     \brief Data structure defining the fields in the nifti2 header.
%            This binary header should be found at the beginning of a valid
%            NIFTI-2 header file.
%  */
%                         /*************************/  /************************/ /************/
% struct nifti_2_header { /* NIFTI-2 usage         */  /* NIFTI-1 usage        */ /*  offset  */
%                         /*************************/  /************************/ /************/
% int   sizeof_hdr;     /*!< MUST be 540           */  /* int sizeof_hdr; (348) */  /*   0 */
% char  magic[8] ;      /*!< MUST be valid signature. */  /* char magic[4];     */  /*   4 */
% short datatype;       /*!< Defines data type!    */  /* short datatype;       */  /*  12 */
% short bitpix;         /*!< Number bits/voxel.    */  /* short bitpix;         */  /*  14 */
% int64_t dim[8];       /*!< Data array dimensions.*/  /* short dim[8];         */  /*  16 */
% double intent_p1 ;    /*!< 1st intent parameter. */  /* float intent_p1;      */  /*  80 */
% double intent_p2 ;    /*!< 2nd intent parameter. */  /* float intent_p2;      */  /*  88 */
% double intent_p3 ;    /*!< 3rd intent parameter. */  /* float intent_p3;      */  /*  96 */
% double pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];      */  /* 104 */
% int64_t vox_offset;   /*!< Offset into .nii file */  /* float vox_offset;     */  /* 168 */
% double scl_slope ;    /*!< Data scaling: slope.  */  /* float scl_slope;      */  /* 176 */
% double scl_inter ;    /*!< Data scaling: offset. */  /* float scl_inter;      */  /* 184 */
% double cal_max;       /*!< Max display intensity */  /* float cal_max;        */  /* 192 */
% double cal_min;       /*!< Min display intensity */  /* float cal_min;        */  /* 200 */
% double slice_duration;/*!< Time for 1 slice.     */  /* float slice_duration; */  /* 208 */
% double toffset;       /*!< Time axis shift.      */  /* float toffset;        */  /* 216 */
% int64_t slice_start;  /*!< First slice index.    */  /* short slice_start;    */  /* 224 */
% int64_t slice_end;    /*!< Last slice index.     */  /* short slice_end;      */  /* 232 */
% char  descrip[80];    /*!< any text you like.    */  /* char descrip[80];     */  /* 240 */
% char  aux_file[24];   /*!< auxiliary filename.   */  /* char aux_file[24];    */  /* 320 */
% int qform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short qform_code;      */  /* 344 */
% int sform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short sform_code;      */  /* 348 */
% double quatern_b ;    /*!< Quaternion b param.   */ /* float quatern_b;       */  /* 352 */
% double quatern_c ;    /*!< Quaternion c param.   */ /* float quatern_c;       */  /* 360 */
% double quatern_d ;    /*!< Quaternion d param.   */ /* float quatern_d;       */  /* 368 */
% double qoffset_x ;    /*!< Quaternion x shift.   */ /* float qoffset_x;       */  /* 376 */
% double qoffset_y ;    /*!< Quaternion y shift.   */ /* float qoffset_y;       */  /* 384 */
% double qoffset_z ;    /*!< Quaternion z shift.   */ /* float qoffset_z;       */  /* 392 */
% double srow_x[4] ;    /*!< 1st row affine transform. */  /* float srow_x[4];  */  /* 400 */
% double srow_y[4] ;    /*!< 2nd row affine transform. */  /* float srow_y[4];  */  /* 432 */
% double srow_z[4] ;    /*!< 3rd row affine transform. */  /* float srow_z[4];  */  /* 464 */
% int slice_code ;      /*!< Slice timing order.   */  /* char slice_code;      */  /* 496 */
% int xyzt_units ;      /*!< Units of pixdim[1..4] */  /* char xyzt_units;      */  /* 500 */
% int intent_code ;     /*!< NIFTI_INTENT_* code.  */  /* short intent_code;    */  /* 504 */
% char intent_name[16]; /*!< 'name' or meaning of data. */ /* char intent_name[16]; */  /* 508 */
% char dim_info;        /*!< MRI slice ordering.   */      /* char dim_info;        */  /* 524 */
% char unused_str[15];  /*!< unused, filled with \0 */                                  /* 525 */
% } ;                   /**** 540 bytes total ****/
% typedef struct nifti_2_header nifti_2_header ;
 
%         sizeof_hdr: 540
%              magic: 'n+2 ←↵'
%           datatype: 32
%             bitpix: 64
%                dim: [6 1 1 1 4096 16 128 1]
%          intent_p1: 0
%          intent_p2: 0
%          intent_p3: 0
%             pixdim: [1 30 40 30 2.5000e-04 1 1 1]
%         vox_offset: 1680
%          scl_slope: 1
%          scl_inter: 0
%            cal_max: 0
%            cal_min: 0
%     slice_duration: 0
%            toffset: 0
%        slice_start: 0
%          slice_end: 0
%            descrip: ''
%           aux_file: ''
%         qform_code: 0
%         sform_code: 2
%          quatern_b: 1
%          quatern_c: 0
%          quatern_d: 0
%          qoffset_x: -1.1622
%          qoffset_y: -11.5254
%          qoffset_z: 16.9492
%             srow_x: [30 0 0 -1.1622]
%             srow_y: [0 -40 0 -11.5254]
%             srow_z: [0 0 -30 16.9492]
%         slice_code: 0
%         xyzt_units: 0
%        intent_code: 0
%        intent_name: 'mrs_v0_5'
%           dim_info: 0
%         unused_str: ''
%          extension: [1 0 0 0]
%            version: 2
%        swap_endian: 0
%          file_name: '/Users/competer/Documents/data_MR_scans/2022-04-20_JNEA1U_HERMESINVIVO/MEGA-GABA/peter_test_spec2nii.nii.gz'

%Default pixdim x,y,z sizes are 10m == no localization, or poor coil
%sensitivity

nii_mrs=struct();

header.sizeof_hdr=540;
%header.magic=magic(8);
header.magic=sprintf('n+%g%s', 2, char([0 13 10 26 10]));
header.datatype=32;
header.bitpix=64;
tmp_dim=4; %Starting tmp_dim(1), minimum dimension is 4 [kx ky kz t] - **PT**2023
if isfield(in.dims,'kx') %MRSI
    if ~in.dims.kx
        kx=1;
    else
        kx=in.sz(in.dims.kx);
    end
    if ~in.dims.ky
        ky=1;
    else
        ky=in.sz(in.dims.ky);
    end
    if ~in.dims.kz
        kz=1;
    else
        kz=in.sz(in.dims.kz);
    end
    dwelltime=in.spectralDwellTime;
else %SVS/FID
    kx=1;
    ky=1;
    kz=1;
    dwelltime=in.dwelltime;
end
tmp_dim(2:4)=[kx,ky,kz];
t=in.sz(in.dims.t);
tmp_dim(5)=t;
dim_cnt=5;

if in.dims.coils
    coils=in.sz(in.dims.coils);
    dim_cnt=dim_cnt+1;
    tmp_dim(dim_cnt)=coils;
    tmp_dim(1)=tmp_dim(1)+1;
end

if in.dims.averages
    averages=in.sz(in.dims.averages);
    dim_cnt=dim_cnt+1;
    tmp_dim(dim_cnt)=averages;
    tmp_dim(1)=tmp_dim(1)+1;
end

if in.dims.subSpecs
    subSpecs=in.sz(in.dims.subSpecs);
    dim_cnt=dim_cnt+1;
    tmp_dim(dim_cnt)=subSpecs;
    tmp_dim(1)=tmp_dim(1)+1;
end

fids=zeros(kx,ky,kz,size(in.fids,in.dims.t));
fids(1,1,1,:)=in.fids;
nii=nii_tool('init',fids);

%dim(1) needs to be +1 from the actual number of dims. Order is: 
%[dim1 kx ky kz t coils averages subspecs extras] - will need to
%incorporate extras dimensions too when I encounter it - **PT**2023
header.dim=tmp_dim;
header.intent_p1=0;
header.intent_p2=0;
header.intent_p3=0;
header.pixdim=[1 in.hdr.MeasYaps.sSpecPara.sVoI.dReadoutFOV ...
    in.hdr.MeasYaps.sSpecPara.sVoI.dPhaseFOV ...
    in.hdr.MeasYaps.sSpecPara.sVoI.dThickness ...
    dwelltime 1 1 1];
header.vox_offset=1680;
header.scl_slope=1; 
header.scl_inter=0;
header.cal_max=0;
header.cal_min=0;
header.slice_duration=0;
header.toffset=0;
header.slice_start=0;
header.slice_end=0;
header.descrip='';
header.aux_file='';
header.qform_code=0;
header.sform_code=2;
header.quatern_b=1;
header.quatern_c=0;
header.quatern_d=0;
header.qoffset_x=-in.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dSag;
header.qoffset_y=-in.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dCor;
header.qoffset_z=in.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dTra;
header.srow_x=[in.hdr.MeasYaps.sSpecPara.sVoI.dReadoutFOV 0 0 header.qoffset_x];
header.srow_y=[0 -in.hdr.MeasYaps.sSpecPara.sVoI.dPhaseFOV 0 header.qoffset_y];
header.srow_z=[0 0 -in.hdr.MeasYaps.sSpecPara.sVoI.dThickness header.qoffset_z];
header.slice_code=0;
header.xyzt_units=0;
header.intent_code=0;
header.intent_name='mrs_v0_5';
header.dim_info=0;
header.unused_str='';
header.extension=[1 0 0 0];
% header.version=2;
% header.swap_endian=0;
% header.filename=[mrs_struct.filename(1:end-2) 'nii.gz'];

header_ext.SpectrometerFrequency=in.txfrq/1e6;
header_ext.ResonantNucleus={in.nucleus};
if in.dims.coils
    header_ext.dim_5='DIM_COIL';
end
if in.dims.averages
    header_ext.dim_6='DIM_DYN';
end
if in.dims.subSpecs
    header_ext.dim_7='DIM_EDIT';
end
header_ext.EchoTime=in.te/1e6;
header_ext.RepetitionTime=in.tr/1e6;
header_ext.ExcitationFlipAngle=in.hdr.Meas.FlipAngle;
header_ext.TxOffset=in.hdr.Meas.dDeltaFrequency;
header_ext.Manufacturer=in.hdr.Meas.Manufacturer;
header_ext.ManufacturersModelName=in.hdr.Meas.ManufacturersModelName;
header_ext.DeviceSerialNumber=in.hdr.Meas.DeviceSerialNumber;
header_ext.SoftwareVersions=in.hdr.Meas.SoftwareVersions;
header_ext.InstitutionName=in.hdr.Meas.InstitutionName;
header_ext.InstitutionAddress=in.hdr.Meas.InstitutionAddress;

if isfield(in.hdr.MeasYaps,'sCoilSelectMeas')
    header_ext.RxCoil=in.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID;
else %different field for older data
    header_ext.RxCoil=in.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID;
end
header_ext.SequenceName=in.hdr.Meas.SequenceString;
header_ext.ProtocolName=in.hdr.Meas.ProtocolName;
header_ext.PatientPosition=in.hdr.Meas.PatientPosition;
header_ext.PatientName='';%mrs_struct.hdr.Meas.PatientName;
header_ext.PatientWeight='';%mrs_struct.hdr.Meas.PatientWeight;
header_ext.PatientDoB='';%mrs_struct.hdr.Meas.PatientBirthday;
header_ext.PatientSex='';%mrs_struct.hdr.Meas.PatientSex;
header_ext.ConversionMethod='FID-A twix2nii v0';
header_ext.ConvertionTime=datetime(now,'ConvertFrom','datenum');
[~,fname,fext]=fileparts(in.filename);
header_ext.OriginalFile=[fname fext];
header_ext.kSpace=[0,0,0];
header_ext.PulseSequenceFile.Value=in.hdr.Meas.tSequenceFileName;
header_ext.PulseSequenceFile.Description='Sequence binary path.';
header_ext.IceProgramFile.Value=in.hdr.Meas.tICEProgramName;
header_ext.IceProgramFile.Description='Reconstruction binary path.';

nii_mrs.hdr=header;
nii_mrs.hdr_ext=header_ext;
out.nii_mrs=nii_mrs;

end