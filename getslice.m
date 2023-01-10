function MRSIStruct=getslice(MRSIStruct)
ori_valid=0;
while ~ori_valid
    user_ori=input('What coordinate? [LR;AP;HF]:    ','s');
    if ~strcmpi(user_ori,{'LR','AP','HF'})
        fprintf('User input not valid [%s]\n',user_ori);
    else
        ori_valid=1;
        switch user_ori
            case {'LR','lr'}
                ori_idx=1;
            case {'AP','ap'}
                ori_idx=2;
            case {'HT','ht'}
                ori_idx=3;
        end
    end
end

user_loc=input('What slice number? [+/-]:');


dcm_path=uigetdir('','Please select scan folder');
dir_dcm_path=dir(dcm_path);

dcm_idx = 3;
tmp_diff = 300;
tmp_full_path = '';
found_img = 0;

while ~found_img
    tmp_info = dicominfo([dcm_path filesep dir_dcm_path(dcm_idx).name]);
    loc_diff = abs(user_loc - tmp_info.ImagePositionPatient(ori_idx));
    if loc_diff<tmp_diff
        tmp_diff = loc_diff;
        dcm_info = tmp_info;
        dcm_idx = dcm_idx +1;
        if dcm_idx>size(dir_dcm_path,1)-2
            disp('Dicom image can not be found');
            return;
        end
    else
        found_img = 1;
        disp(['The matching dicom image is: ' dcm_info.Filename]);
        sprintf('The closest matching slice to (%f [%s]) is: %s',user_loc,user_ori,num2str(tmp_info.ImagePositionPatient(ori_idx))); 
        figure;
        imshow(dicomread(dcm_info),[]);
        MRSIStruct.mrsi.dcm_overlay=dcm_info;
        MRSIStruct.flags.dcmfound=1;
%         dcm_file=dcm_info.Filename;
%         save_dir = uigetdir('/ResearchData4/ptruong/CAMH-DATA/Phosphorus','Select the save directory');
%         unix(['cp ' dcm_info.Filename ' ' save_dir filesep]);
    end
end