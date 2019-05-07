close all; clear all;

%forward solution from OpenMEEG
file_fwd={
    's031_d10-lh.gain.mat';
    's031_d10-rh.gain.mat';
    };

%SEEG contact info
file_seeg_contact_postMR='electrode_040119_085342.mat';

%source space
file_source_fif='/Users/fhlin_admin/workspace/seeg/subjects/s031/bem/s031-5-src.fif';

file_output='seeg_fwd_050719.mat';

%%%%%%%%%%%%%%%%%%%%%%%
[src] = mne_read_source_spaces(file_source_fif);

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
end;

save(file_output,'A');

