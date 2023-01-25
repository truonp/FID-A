function MRSIStruct=op_CSIoverlayMRIslice_LCModel(MRSIStruct)

if ~isfield(MRSIStruct.flags,'dcmfound')
    disp('You must run the function <getslice> first, to get the matching dicom image');
    return
end

sd_thresh=25;%cutoff threshold for %SD

%There will be cases where the CSI FOV is smaller than the image, but focus
%on case where it's larger for now

dcm_img=dicomread(MRSIStruct.mrsi.dcm_overlay);
% [dcm_file dcm_path]=uigetfile('*.IMA');
% dcm_info=dicominfo([dcm_path filesep dcm_file]);
% dcm_info=dicominfo('/Users/competer/Documents/data_MR_scans/2022-09-30_JNEA1U_31P_test/JNEA1U_31P_TEST_JNEA1U_31P_TEST/TECHNICAL_JNEA1U_20220930_121308_959000/COR_0004/JNEA1U_31P_TEST.MR.TECHNICAL_JNEA1U.0004.0048.2022.09.30.13.39.21.928538.15245510.IMA');
% dcm_img=dicomread(dcm_info);

% csi_pixel_xy=[((-csi_av.fov.x/2):csi_av.fov.x/2).*dcm_info.PixelSpacing(1);((-csi_av.fov.y/2):csi_av.fov.y/2).*dcm_info.PixelSpacing(2)]; %csi fov (w.r.t. dicom image mm)
% csi_edges_xy=[csi_pixel_xy(1,1:csi_av.voxelSize.x:end);csi_pixel_xy(2,1:csi_av.voxelSize.y:end)]; %csi grid edges (w.r.t. dicom image mm)
% csi_pixel_xy=[((-csi_av.fov.x/2):dcm_info.PixelSpacing(1):csi_av.fov.x/2);((-csi_av.fov.y/2):dcm_info.PixelSpacing(2):csi_av.fov.y/2)];
% csi_edges_xy=[csi_pixel_xy(1,1:csi_av.voxelSize.x:end);csi_pixel_xy(2,1:csi_av.voxelSize.y:cdend)];
csi_pixel_xy=[((-MRSIStruct.fov.x/2):MRSIStruct.fov.x/2);((-MRSIStruct.fov.y/2):MRSIStruct.fov.y/2)];
csi_edges_xy=[csi_pixel_xy(1,1:MRSIStruct.voxelSize.x:end);csi_pixel_xy(2,1:MRSIStruct.voxelSize.y:end)];

% dcm_fov=[size(dcm_img).*dcm_info.PixelSpacing']; %field of view for the dicom image (mm)
% dcm_pixel_xy=[((-size(dcm_img,1)/2):size(dcm_img,1)/2).*dcm_info.PixelSpacing(1);((-size(dcm_img,2)/2):size(dcm_img,2)/2).*dcm_info.PixelSpacing(2)]; %pixel locations (in mm) of the image
dcm_fov=[size(dcm_img)]; %field of view for the dicom image (mm)
dcm_pixel_xy=[((-size(dcm_img,1)/2):size(dcm_img,1)/2);((-size(dcm_img,2)/2):size(dcm_img,2)/2)]; %pixel locations (in mm) of the image

csi_centre=MRSIStruct.imageOrigin; %(-+)rl;(-+)ap;(-+)fh(is)

tmp_csi_shift=csi_pixel_xy+[-csi_centre(2); csi_centre(3)];
[~,shifted_origin]=min(abs(tmp_csi_shift),[],2);

new_img=zeros(MRSIStruct.fov.x,MRSIStruct.fov.y);
new_img(shifted_origin(1)-size(dcm_img,1)/2+1:shifted_origin(1)+size(dcm_img,1)/2,shifted_origin(2)-size(dcm_img,2)/2+1:shifted_origin(2)+size(dcm_img,2)/2)=dcm_img;
imshow(new_img,[]); hold on; xline(1:MRSIStruct.voxelSize.x:MRSIStruct.fov.x,'g'); yline(1:MRSIStruct.voxelSize.x:MRSIStruct.fov.x,'g');

if ~isfield(MRSIStruct.mrsi,'lcm')
    [csv_name,csv_path]=uigetfile('*.csv');
%     csv_path = '/Users/competer/Documents/data_MR_scans/2022-09-30_JNEA1U_31P_test/FIDA_CSI_output';
%     csv_name = 'test.csv';
    csv_fullpath=fullfile(csv_path,csv_name);
    metab_data = readtable(csv_fullpath);
    MRSIStruct.mrsi.lcm.csv_file=csv_fullpath;
    MRSIStruct.mrsi.lcm.met_data=metab_data;
else
    metab_data=MRSIStruct.mrsi.lcm.met_data;
end

lcm_names = metab_data.Properties.VariableNames(3:end);

y_range = [min(metab_data.Row),max(metab_data.Row)];
x_range = [min(metab_data.Col),max(metab_data.Col)];
metab_arr = zeros([size(metab_data,1),1]); %metabolite values
tmp_metab_arr = metab_arr;
metab_mat = zeros(x_range(2)-x_range(1)+1,y_range(2)-y_range(1)+1); %matrix grid mapping met_arr values to x,y coordinates

sd_arr1=zeros(size(metab_data,1));
sd_arr2=sd_arr1;

n_metab=0;
tmp_str='';
for ii=1:3:size(lcm_names,2)
    n_metab=n_metab+1;
    tmp_str=sprintf('%s[%d]: %s\n',tmp_str,n_metab,lcm_names{ii});
end
metab_list=tmp_str(1:end-1);

fprintf('Which metabolite do you wish to display?');

u_valid=0;
while ~u_valid
    metab_num=input(sprintf('\n%s\n',metab_list));
    if metab_num < 1 || metab_num > n_metab
        fprintf('Input is out of bounds [1,%d]',n_metab);
    else
        u_valid=1;
        sd_arr1=metab_data.(lcm_names{metab_num*3-1});
        tmp_metab_arr=metab_data.(lcm_names{metab_num*3-2}).*(~any( sd_arr1 >sd_thresh,2));
    end
end

no_denom=1;
r_valid=0;
while ~r_valid
    include_denom=input(sprintf('Do you wish to show ratios? [Y/N]'),'s');
    if strcmpi(include_denom,'y')
        fprintf('Which metabolite do you wish to have as denominator?');
        u_valid=0;
        while ~u_valid
            metab_denom=input(sprintf('\n%s\n',metab_list));
            if metab_denom < 1 || metab_denom > n_metab
                fprintf('Input is out of bounds [1,%d]',n_metab);
            else
                u_valid=1;
                r_valid=1;
                no_denom=0;
                sd_arr2=metab_data.(lcm_names{metab_denom*3-1});
                tmp_metab_arr=tmp_metab_arr.*(~any( sd_arr2 >sd_thresh,2))./metab_data.(lcm_names{3*metab_denom-2});   
                tmp_metab_arr(isnan(tmp_metab_arr))=0;
                metab_arr=tmp_metab_arr.*(~any( isinf( tmp_metab_arr),2));
            end
        end
    elseif strcmpi(include_denom,'n')
        r_valid=1;
        metab_arr=tmp_metab_arr.*(~any( isinf( tmp_metab_arr) | isnan( tmp_metab_arr),2));
    else
        sprintf('Invalid response, please input either [Y/N]');
    end
end
    
for i = 1:length(metab_arr)
    x = metab_data.Col(i);
    y = metab_data.Row(i);
    metab_mat(x,y) = metab_arr(i);
end

met_img = imresize(metab_mat,size(new_img),'nearest');

%snippet of this code taken from: https://www.mathworks.com/matlabcentral/answers/635094-overlay-the-image-with-transparency
cmap = turbo(256);
rgbImage = ind2rgb(uint8(255 * mat2gray(met_img)), cmap);
ha = imshow(new_img, []);
hold on;
hb = imshow(rgbImage);hb.AlphaData = 0.3;
end

