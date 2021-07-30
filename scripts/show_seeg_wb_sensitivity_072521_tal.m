close all; clear all;
clear global  etc_render_fsbrain;

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects/'); 
subject='fsaverage';
surf='orig';

% mri=MRIread('/Users/fhlin/workspace/seeg/subjects/s036/mri/orig.mgz'); %for MAC/Linux
mri=MRIread('/Users/fhlin/workspace/seeg/subjects/fsaverage/mri/orig.mgz'); %for MAC/Linux
% 
% %volume stc
% 
mri_overlay=MRIread('average_seeg_wb_sensitiivty_072521_snr-vol-tal.2mm.mgh');

targ=MRIread('/Applications/freesurfer/average/mni305.cor.subfov2.mgz'); %MNI-Talairach space with 2mm resolution (for MAC)
%targ=MRIread(sprintf('%s/average/mni305.cor.subfov2.mgz',getenv('FREESURFER_HOME'))); %MNI-Talairach space with 2mm resolution (for server)

targ_reg=etc_read_xfm('file_xfm','/Applications/freesurfer/average/mni305.cor.subfov2.reg'); %MNI-Talairach space with 2mm resolution (for MAC)
%targ_reg=etc_read_xfm('file_xfm',sprintf('%s/average/mni305.cor.subfov2.reg',getenv('FREESURFER_HOME'))); %MNI-Talairach space with 2mm resolution (for server)

etc_render_fsbrain('surf','inflated','hemi','lh','subject',subject,'vol_reg',targ_reg,'vol',targ,'overlay_vol',mri_overlay,'overlay_threshold',[0.01 1].*1e-10);

