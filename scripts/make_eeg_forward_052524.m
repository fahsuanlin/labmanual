close all; clear all;

%forward solution from OpenMEEG
file_fwd={
    'eeg_lh_gain.mat';
    'eeg_rh_gain.mat';
    };

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri/subjects/'); %for MAC/Linux

subject='180411_PYY';

file_eeg_mat='./fwd_prep_041618.mat';

%source space
file_source_fif=sprintf('%s/%s/bem/%s-5-src.fif',getenv('SUBJECTS_DIR'),subject,subject);
select_channel={ 'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'};

file_output='eeg_fwd_052524.mat';

%%%
[src] = mne_read_source_spaces(file_source_fif);


load(file_eeg_mat); %including variable 'n_meg'
[A_label,i1,i2]=intersect(lower(select_channel),lower(points_label));
 

for hemi_idx=1:2
    load(file_fwd{hemi_idx});

    A(hemi_idx).label=A_label;
    for v_idx=1:length(A_label)
        A(hemi_idx).coord(v_idx,:)=points(i2(v_idx),1:3).*1e3;
    end;

    A(hemi_idx).A=linop;

    A(hemi_idx).v_idx=double(src(hemi_idx).vertno-1);
    A(hemi_idx).ori=double(src(hemi_idx).nn(src(hemi_idx).vertno,:));
    A(hemi_idx).loc=double(src(hemi_idx).rr(src(hemi_idx).vertno,:));
    A(hemi_idx).source_file=file_source_fif;
end;

save(file_output,'A');
