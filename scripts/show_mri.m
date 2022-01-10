close all; clear all;

subject='s031';
hemi='rh';
surf='inflated';


file_overlay_register='../fmri_data/unpack/register.dat';
file_overlay_vol='../fmri_data/unpack/bold/006/f.mgz';
overlay_threshold=[100 1500];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects/'); %for MAC/Linux
%setenv('SUBJECTS_DIR','D:\fhlin\Users\fhlin\workspace\seeg\subjects'); %for PC

mri=MRIread(sprintf('/Users/fhlin/workspace/seeg/subjects/%s/mri/orig.mgz',subject)); %for MAC/Linux
%mri=etc_MRIread('D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\orig.mgz'); %for PC

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm',sprintf('%s/%s/mri/transforms/talairach.xfm',getenv('SUBJECTS_DIR'),subject)); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\transforms\talairach.xfm'); %for PC

%%% volume overlay
f=MRIread(file_overlay_vol);
r=etc_read_xfm('file_xfm',file_overlay_register);

%prepare overlay volume matched to the underlay volume
vol = etc_MRIvol2vol(f,mri,r,'frames',[1]); %only the first volume
 
%prepare overlay surface
vol2surf = etc_MRIvol2surf(f,surf,r,'subject',subject,'frames',[1],'hemi',hemi);

etc_render_fsbrain('surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'overlay_vol',vol,'overlay_value',vol2surf,'overlay_vertex',[1:size(vol2surf,1)]-1,'overlay_threshold',overlay_threshold);
