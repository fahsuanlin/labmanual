close all; clear all;
clear global  etc_render_fsbrain;

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/seeg/subjects/'); 
subject='s031';
surf='orig';

file_forward_mat='seeg_fwd_wb_091019.mat';
%file_forward_mat='seeg_dec_fwd_wb_091019.mat';

mri=MRIread('/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/orig.mgz'); %for MAC/Linux

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/transforms/talairach.xfm'); %for MAC/Linux

load(file_forward_mat);
%append the forward matrix object 'A' surface coordinates; used for
%interpolation
hemi={'lh','rh'};
for hemi_idx=1:2
   file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hemi{hemi_idx},surf);
   [A(hemi_idx).vertex_coords, A(hemi_idx).faces] = read_surf(file_surf);
end;

%volume stc
[vol_stc,vv,d0,d1,timeVec]=inverse_read_stc('seeg_wb_mne_091019_v-vol.stc');


%without specified electrode
etc_render_fsbrain('surf','orig','hemi','rh','subject','s031','vol',mri,'talxfm',(talxfm),'overlay_vol_stc',vol_stc,'vol_A',A,'overlay_stc_timeVec',timeVec);

