close all; clear all;

targ=MRIread('/Users/fhlin_admin/workspace/7t_music_skku/subjects/SUB3B/mri/orig.mgz');
subject='SUB3B';

targ_reg=eye(4);

dti=MRIread('../unpack/bold/010/fmcprstc.nii.gz');


% do registration between pre- and post-OP by the following command:
%
% cd /Users/fhlin_admin/workspace/seeg/subjects/s057/tmp
% bbregister --s s057 --mov ../../../s057/DTI/nodif_brain.nii.gz --init-coreg --reg register.dat --t2
% 
% %check registration
% tkregisterfv --mov ../../../s057/DTI/nodif_brain.nii.gz --reg register.dat --surfs

%dti_reg=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/7t_music_skku/subjects/SUB4A/tmp/register.dat'); %for MAC/Linux
%dti_reg=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\s057\tmp\register.dat'); %for PC

dti_reg=etc_read_xfm('file_xfm','bb_register_init.dat'); %read the pre-defined registration matrix

dtim=MRIvol2vol(dti,targ,dti_reg); %transform the image to match the structural MRI volume

%%%%% if needed, use the registration toool ('k' in etc_render_fsbrain) to
%%%%% register manually and expert the registration matrix as
%%%%% `overlay_xfm`. Then use the product dti_reg*overlay_xfm as the
%%%%% transformation
% dtim=MRIvol2vol(dti,targ,dti_reg*overaly_xfm);

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/7t_music_skku/subjects/SUB4A/mri/transforms/talairach.xfm'); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\s057\mri\transforms\talairach.xfm'); %for PC

etc_render_fsbrain('surf','orig','hemi','rh','subject',subject,'vol',targ,'overlay_vol',dtim,'overlay_threshold',[100 500]); 

return;
