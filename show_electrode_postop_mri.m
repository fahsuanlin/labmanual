close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/seeg/subjects/'); %for MAC/Linux
%setenv('SUBJECTS_DIR','D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects'); %for PC

mri_post=MRIread('/Users/fhlin_admin/workspace/seeg/subjects/s031_post/mri/orig.mgz'); %for MAC/Linux
%mri_post=etc_MRIread('D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\2036_post\mri\orig.mgz'); %for PC

mri=MRIread('/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/orig.mgz'); %for MAC/Linux
%mri=etc_MRIread('D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\2036\mri\orig.mgz'); %for PC

% do registration between pre- and post-OP by the following command:
%
% cd /Users/fhlin_admin/workspace/seeg/subjects/s031_post/tmp
% bbregister --s s031 --mov ../mri/orig.mgz --init-fsl --reg register.dat --t1
%
xfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/seeg/subjects/s031_post/tmp/register.dat'); %for MAC/Linux
%xfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\2036_post\tmp\register.dat'); %for PC

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/transforms/talairach.xfm'); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin_admin\workspace\seeg\subjects\2036\mri\transforms\talairach.xfm'); %for PC

etc_render_fsbrain('surf','inflated','hemi','rh','subject','s031_post','vol',mri_post,'vol_pre_xfm',inv(xfm)*mri.vox2ras*inv(mri_post.vox2ras),'talxfm',(talxfm));
view(90,30);