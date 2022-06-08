close all; clear all;

%file_stc='erp_mne_outside_060722_003_dspm_cc-lh.stc';
file_stc='erp_mne_outside_060722_003_dspm_cc-rh.stc';

[stc,v,t0,t1,timeVec]=inverse_read_stc(file_stc);

%etc_render_fsbrain('hemi','lh','subject','180322_PYW','overlay_vertex',v,'overlay_stc',stc,'overlay_threshold',[5 15],'overlay_stc_timeVec',timeVec)
etc_render_fsbrain('hemi','rh','subject','180322_PYW','overlay_vertex',v,'overlay_stc',stc,'overlay_threshold',[5 15],'overlay_stc_timeVec',timeVec)

return;

