close all; clear all;

subject='011624';
hemi='rh';
surf='white';


file_overlay_register='./bb_register_blipdown_init.dat';
%file_overlay_vol='epi__epi_blipdown.mgh';
file_overlay_vol='epi__epi_blipdown_fliplr.mgh'; %use this one because SMS-InI has flipped A-P direction compared to EPI.
overlay_threshold=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setenv('SUBJECTS_DIR','/Users/fhlin/workspace/smsini_blipud/subjects/'); %for MAC/Linux
%setenv('SUBJECTS_DIR','D:\fhlin\Users\fhlin\workspace\seeg\subjects'); %for PC

mri=MRIread(sprintf('/Users/fhlin/workspace/smsini_blipud/subjects/%s/mri/orig.mgz',subject)); %for MAC/Linux
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

%etc_render_fsbrain('surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'overlay_vol',vol,'overlay_value',vol2surf,'overlay_vertex',[1:size(vol2surf,1)]-1,'overlay_threshold',overlay_threshold,'overlay_flag_paint_on_cortex',0);
etc_render_fsbrain('surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'overlay_vol',vol,'overlay_threshold',overlay_threshold,'overlay_flag_paint_on_cortex',0);

%1. calcualte r * overlay_xfm after exporting the registration in Matlab work space
%
%2. check initial registration at command line: tkregisterfv --mov smsini_mb_run_1_ref.mgh --reg bb_register_init.dat --surfs  --sd /Users/fhlin/workspace/smsini_blipud/subjects/
