close all; clear all;

%forward solution from OpenMEEG
file_fwd={
    'seeg027_d10_wb_091019-lh.gain.mat';
    'seeg027_d10_wb_091019-rh.gain.mat';
    };

%SEEG contact info
file_seeg_contact_postMR='electrode_081420_152935.mat';

%source space
%file_source_fif='/space_lin2/fhlin/seeg/subjects/seeg027/bem/seeg027-5-src.fif';
%file_source_wholebrain='/space_lin2/fhlin/seeg/subjects/seeg027/mri/aseg.mgz';
%file_source_wholebrain_orig='/space_lin2/fhlin/seeg/subjects/seeg027/mri/orig.mgz';


file_source_fif='/space_lin2/fhlin/seeg_ictal_source/subjects/seeg027/bem/seeg027-5-src.fif';
file_source_wholebrain='/space_lin2/fhlin/seeg_ictal_source/subjects/seeg027/mri/aseg.mgz';
file_source_wholebrain_orig='/space_lin2/fhlin/seeg_ictal_source/subjects/seeg027/mri/orig.mgz';

wholebrain_index={
    [9 10 11 12 13 16 17 18 19 20 26 27 ];
    [48 49 50 51 52 53 54 55 56 58 59];
    };
%the following indices were from FreeSurfer look-up table 
% 9   Left-Thalamus                           0   118 14  0
% 10  Left-Thalamus-Proper*                   0   118 14  0
% 11  Left-Caudate                            122 186 220 0
% 12  Left-Putamen                            236 13  176 0
% 13  Left-Pallidum                           12  48  255 0
% 16  Brain-Stem                              119 159 176 0
% 17  Left-Hippocampus                        220 216 20  0
% 18  Left-Amygdala                           103 255 255 0
% 19  Left-Insula                             80  196 98  0
% 20  Left-Operculum                          60  58  210 0
% 26  Left-Accumbens-area                     255 165 0   0
% 27  Left-Substancia-Nigra                   0   255 127 0
% 48  Right-Thalamus                          0   118 14  0
% 49  Right-Thalamus-Proper*                  0   118 14  0
% 50  Right-Caudate                           122 186 220 0
% 51  Right-Putamen                           236 13  176 0
% 52  Right-Pallidum                          13  48  255 0
% 53  Right-Hippocampus                       220 216 20  0
% 54  Right-Amygdala                          103 255 255 0
% 55  Right-Insula                            80  196 98  0
% 56  Right-Operculum                         60  58  210 0
% 58  Right-Accumbens-area                    255 165 0   0
% 59  Right-Substancia-Nigra                  0   255 127 0

file_output='seeg_fwd_wb_091019.mat';

%%%%%%%%%%%%%%%%%%%%%%%
[src] = mne_read_source_spaces(file_source_fif);
%whole brain source space
if(~isempty(file_source_wholebrain_orig))
    src_wb_orig=MRIread(file_source_wholebrain_orig);
end;
if(~isempty(file_source_wholebrain))
    src_wb=MRIread(file_source_wholebrain);
    src_wb_vol=src_wb.vol;
    
    for hemi_idx=1:2
        src_wb_idx{hemi_idx}=[];
        wb_index{hemi_idx}=[];
        for idx=1:length(wholebrain_index{hemi_idx})
            tmp=find(src_wb_vol(:)==wholebrain_index{hemi_idx}(idx));
            src_wb_idx{hemi_idx}=cat(1,src_wb_idx{hemi_idx},tmp(:));
            wb_index{hemi_idx}=cat(1,wb_index{hemi_idx},ones(size(tmp(:))).*wholebrain_index{hemi_idx}(idx));
        end;
        
        for idx=1:length(src_wb_idx{hemi_idx})
            [cc,rr,ss]=ind2sub(size(src_wb_orig.vol),src_wb_idx{hemi_idx}(idx));
            tmp=src_wb.tkrvox2ras*[rr cc ss 1]';
            wb_coord{hemi_idx}(idx,:)=tmp(1:3)'; %this is the coordinate of whole-brain source in the surface coordinates
        end;

        [cc,rr,ss]=ind2sub(size(src_wb_orig.vol),src_wb_idx{hemi_idx}(:));
        cc_min=min(cc(:));
        cc_max=max(cc(:));
        rr_min=min(rr(:));
        rr_max=max(rr(:));
        ss_min=min(ss(:));
        ss_max=max(ss(:));
        cc_idx=ismember(cc(:),[cc_min:2:cc_max]);
        rr_idx=ismember(rr(:),[rr_min:2:rr_max]);
        ss_idx=ismember(ss(:),[ss_min:2:ss_max]);
        crs_idx=(cc_idx&rr_idx&ss_idx);
        src_wb_idx{hemi_idx}=src_wb_idx{hemi_idx}(crs_idx);

        wb_coord_full{hemi_idx}=wb_coord{hemi_idx};
        wb_index_full{hemi_idx}=wb_index{hemi_idx};

        wb_coord{hemi_idx}=wb_coord{hemi_idx}(crs_idx,:);
        wb_index{hemi_idx}=wb_index{hemi_idx}(crs_idx);
    end;
end;


load(file_seeg_contact_postMR);
name={};
for e_idx=1:length(electrode)
    for c_idx=1:electrode(e_idx).n_contact
        name{end+1}=sprintf('%s%d',electrode(e_idx).name,c_idx);
    end;
end;

for hemi_idx=1:2
    load(file_fwd{hemi_idx});
       
    A(hemi_idx).A=linop;
    
    A(hemi_idx).v_idx=double(src(hemi_idx).vertno-1);
    A(hemi_idx).ori=double(src(hemi_idx).nn(src(hemi_idx).vertno,:));
    A(hemi_idx).loc=double(src(hemi_idx).rr(src(hemi_idx).vertno,:));
    A(hemi_idx).source_file=file_source_fif;
    A(hemi_idx).name=name;
    
    A(hemi_idx).src_wb_idx=src_wb_idx{hemi_idx};
    A(hemi_idx).wb_loc_full=wb_coord_full{hemi_idx}./1e3; %coordinates
    A(hemi_idx).wb_index_full=wb_index_full{hemi_idx}; %label based on FreeSurfer LUT
    
    A(hemi_idx).wb_loc=wb_coord{hemi_idx}./1e3; %coordinates
    A(hemi_idx).wb_index=wb_index{hemi_idx}; %label based on FreeSurfer LUT
    A(hemi_idx).wb_source_file=file_source_wholebrain;
end;

save(file_output,'A');

