function MRSIStruct=calculate_ph(MRSIStruct)

coord_path=uigetdir;
dir_coord=dir([coord_path filesep '*.coord']);

for ii=1:size(dir_coord,1)
    %     [coord_name coord_path]=uigetfile('*.coord');
    %     coord_name='test_sl1_6-6.coord';
    %     coord_path='/Users/competer/Documents/data_MR_scans/2022-09-30_JNEA1U_31P_test/FIDA_CSI_output';
    coord_fullpath=fullfile(dir_coord(ii).folder,dir_coord(ii).name);
    raw_coord=readlines(coord_fullpath);
    
    coord_file_mat=strings(MRSIStruct.sz(2),MRSIStruct.sz(3));
    pH_mat=zeros(MRSIStruct.sz(2),MRSIStruct.sz(3));
    if ~contains(raw_coord,'FATAL')
        
        [~,coord_name,~]=fileparts(coord_fullpath);
        underscore_find=strfind(coord_name,'_');
        coord_str=coord_name(underscore_find(end)+1:end);
        
        dash_find=strfind(coord_str,'-');
        row=coord_str(1:dash_find-1);
        col=coord_str(dash_find+1:end);
        coord_field_str=sprintf('r%s_r%s',[repmat('0',2-size(row)) row],[repmat('0',2-size(col)) col]);
        
        
        index_start.ppm=find(contains(raw_coord,'points on ppm-axis'));
        index_start.spec=find(contains(raw_coord,'NY phased data points follow'));
        
        tmp_idx=0;
        tmp_idx=find(contains(raw_coord,'PCr'));
        index_start.PCr=tmp_idx(2);
        PCr_conc=raw_coord{tmp_idx(1)};
        perc_find=strfind(PCr_conc,'%');
        
        %getting the SD%, if it is higher than 20%, then don't use data
        %point
        if (str2num(PCr_conc(perc_find-3:perc_find-1))>20)
            fprintf('%s - Poor fit for PCr: %d%%\n',coord_field_str,str2num(PCr_conc(perc_find-3:perc_find-1)));
            continue;
        end
        
        tmp_idx=find(contains(raw_coord,'ATPa'));
        index_start.ATPa=tmp_idx(2);
        
        tmp_idx=find(contains(raw_coord,'Pi'));
        index_start.Pi=tmp_idx(2);
        Pi_conc=raw_coord{tmp_idx(1)};
        %getting the SD%, if it is higher than 25%, then don't use data
        %point
        if (str2num(Pi_conc(perc_find-3:perc_find-1))>25)
            fprintf('%s Poor fit for Pi: %d%%\n',coord_field_str,str2num(Pi_conc(perc_find-3:perc_find-1)));
            continue;
        end
        
        tmp_idx=find(contains(raw_coord,'NADH'));
        index_start.NADH=tmp_idx(2);
        
        ppm_arr=str2num(strjoin(raw_coord(index_start.ppm+1:index_start.spec-1)));
        PCr_arr=str2num(strjoin(raw_coord(index_start.PCr+1:index_start.ATPa-1)));
        Pi_arr=str2num(strjoin(raw_coord(index_start.Pi+1:index_start.NADH-1)));
        
        [~,PCr_idx]=max(PCr_arr);
        [~,Pi_idx]=max(Pi_arr);
        
        dppm=ppm_arr(Pi_idx)-ppm_arr(PCr_idx);
        pKa=6.77;
        dHA=3.23;
        dA=5.7;
        pH=pKa+log((dppm-dHA)/(dA-dppm));
        
        pH_mat(row,col)=pH;
        coord_file_mat(row,col)=coord_fullpath;

    end
end
end