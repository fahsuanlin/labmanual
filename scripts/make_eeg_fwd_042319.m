close all; clear all;

%forward solution from OpenMEEG
file_fwd={
    'eeg_lh_gain.mat';
    'eeg_rh_gain.mat';
    };

%source space
file_source_fif='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/041019-5-src.fif';

file_output='eeg_fwd_042319.mat';

%%%
[src] = mne_read_source_spaces(file_source_fif);

for hemi_idx=1:2
    load(file_fwd{hemi_idx});
       
    A(hemi_idx).A=linop;
    
    A(hemi_idx).v_idx=double(src(hemi_idx).vertno-1);
    A(hemi_idx).ori=double(src(hemi_idx).nn(src(hemi_idx).vertno,:));
    A(hemi_idx).loc=double(src(hemi_idx).rr(src(hemi_idx).vertno,:));
    A(hemi_idx).source_file=file_source_fif;
end;

save(file_output,'A');

