close all; clear all;

%forward solution from OpenMEEG
file_fwd={
    '180322_PYW_d10_wb_eeg_060722-lh.gain.mat';
    '180322_PYW_d10_wb_eeg_060722-rh.gain.mat';
    };

file_source_fif='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/180322_PYW-5-src.fif';
file_source_wholebrain='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/mri/aseg.mgz';
file_source_wholebrain_orig='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/mri/orig.mgz';

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

%sesnor info
file_eeg_mat='./eeg_fwd_prep_041721.mat';
select_channel={ 'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'};

file_output='eeg_fwd_wb_060722.mat';

%%%%%%%%%%%%%%%%%%%%%%%
load(file_eeg_mat);

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
        wb_coord{hemi_idx}=wb_coord{hemi_idx}(crs_idx,:);
%         rr=rr(find(crs_idx));
%         cc=cc(find(crs_idx));
%         ss=ss(find(crs_idx));

    end;
end;


Index=find(contains(points_label,select_channel));
if(length(Index)<=length(select_channel)) %all electrodes were found on topology
    for ii=1:length(select_channel)
        for idx=1:length(points_label)
            if(strcmp(points_label{idx},select_channel{ii}))
                Index(ii)=idx;
                electrode_data_idx(idx)=ii;
            end;
        end;
    end;
    if(isempty(Index))
        fprintf('cannot find corresponding channels!\nerror in loading the topology!\n');
        
        return;
    end;
end;
                 

for hemi_idx=1:2
    load(file_fwd{hemi_idx});
       
    A(hemi_idx).A=linop;
    
    A(hemi_idx).v_idx=double(src(hemi_idx).vertno-1);
    A(hemi_idx).ori=double(src(hemi_idx).nn(src(hemi_idx).vertno,:));
    A(hemi_idx).loc=double(src(hemi_idx).rr(src(hemi_idx).vertno,:));
    A(hemi_idx).source_file=file_source_fif;
    
    A(hemi_idx).src_wb_idx=src_wb_idx{hemi_idx};
    for ch_idx=1:n_eeg
        A(hemi_idx).name{ch_idx}=points_label{Index(ch_idx)};
    end;

    A(hemi_idx).wb_loc=wb_coord{hemi_idx}./1e3; %coordinates
    A(hemi_idx).wb_index=wb_index{hemi_idx}; %label based on FreeSurfer LUT
    A(hemi_idx).wb_source_file=file_source_wholebrain;

end;





save(file_output,'A');