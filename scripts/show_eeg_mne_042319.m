close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/sinica_meg/subjects/'); %for MAC/Linux

subject='041019';
%hemi='lh';
hemi='rh';
threshold=[5 30];

%[stc,v_idx,a,b,timeVec]=inverse_read_stc('eeg_042319_120-lh.stc');
[stc,v_idx,a,b,timeVec]=inverse_read_stc('eeg_042319_120-rh.stc');

etc_render_fsbrain('subject',subject,'hemi',hemi,'overlay_stc',stc,'overlay_vertex',v_idx,'overlay_threshold',threshold,'overlay_smooth',3,'overlay_stc_timeVec',timeVec);
