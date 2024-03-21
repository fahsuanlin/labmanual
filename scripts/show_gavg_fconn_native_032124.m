close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects');

subject='s026';

f=MRIread('hippo_fconn_native_vol_aseg_032124_gavg_hippo_left-anat.nii');
etc_render_fsbrain('subject',subject,'overlay_vol',f,'overlay_threshold',[2 5],'overlay_truncate_neg',1);
