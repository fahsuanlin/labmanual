close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/seeg/subjects/'); %for MAC/Linux

subject='s031';
%hemi='lh';
hemi='rh';
threshold=[30 50];

%[stc,v_idx,a,b,timeVec]=inverse_read_stc('seeg_mne_050719_a-lh.stc');
[stc,v_idx,a,b,timeVec]=inverse_read_stc('seeg_mne_050719_a-rh.stc');

etc_render_fsbrain('subject',subject,'hemi',hemi,'overlay_stc',stc,'overlay_vertex',v_idx,'overlay_threshold',threshold,'overlay_smooth',5,'overlay_stc_timeVec',timeVec);
