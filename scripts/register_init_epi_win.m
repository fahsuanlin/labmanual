close all; clear all;

subject='s002';
hemi='rh';
surf='white';


file_overlay_register='./bb_register_init.dat';
file_overlay_vol='epi_005_f.mgh';
overlay_threshold=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri_music/subjects/'); %for MAC/Linux
setenv('SUBJECTS_DIR','C:\Users\sr810\workspace\eegmri_music\subjects'); %for PC

%mri=MRIread(sprintf('/Users/fhlin/workspace/eegmri_music/subjects/%s/mri/orig.mgh',subject)); %for MAC/Linux
mri=etc_MRIread('C:\Users\sr810\workspace\eegmri_music\subjects\s002\mri\orig.mgh'); %for PC

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

%etc_render_fsbrain('surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'overlay_vol',vol,'overlay_value',vol2surf,'overlay_vertex',[1:size(vol2surf,1)]-1,'overlay_threshold',overlay_threshold,'overlay_flag_paint_on_cortex',0);
etc_render_fsbrain('surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'overlay_vol',vol,'overlay_threshold',overlay_threshold,'overlay_flag_paint_on_cortex',0);
