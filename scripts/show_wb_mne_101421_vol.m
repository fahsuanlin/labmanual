close all; clear all;

file_vol_stc='seeg_wb_mne_cc_raw_101421_lf-vol.stc';

file_forward_mat='seeg_fwd_wb_dec_072521.mat';

subject='s031';

timeVec_range=[-10 500];
hemi='rh';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vol=MRIread('/Users/fhlin/workspace/seeg/subjects/s031/mri/orig.mgz');
vol_post=MRIread('/Users/fhlin/workspace/seeg/subjects/s031_post/mri/orig.mgz');

% do registration between pre- and post-OP by the following command:
%
% cd /Users/fhlin/workspace/seeg/subjects/s031_post/tmp
% bbregister --s s031 --mov ../mri/orig.mgz --init-fsl --reg register.dat --t1
%
xfm=etc_read_xfm('file_xfm','/Users/fhlin/workspace/seeg/subjects/s031_post/tmp/register.dat'); %for MAC/Linux
%xfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\s031_post\tmp\register.dat'); %for PC

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm','/Users/fhlin/workspace/seeg/subjects/s031/mri/transforms/talairach.xfm'); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\s031\mri\transforms\talairach.xfm'); %for PC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert electrode coordinates from post-op MRI to pre-op MRI

file_mat='electrode_040119_085342.mat'; %electrode coordinates in post-op MRI
load(file_mat);

electrode_out=electrode;
for e_idx=1:length(electrode)
    for c_idx=1:electrode(e_idx).n_contact
        surface_coord=electrode(e_idx).coord(c_idx,:);
        surface_coord=inv(xfm)*[surface_coord(:); 1];
        electrode_out(e_idx).coord(c_idx,:)=surface_coord(1:3);
    end;
end;

electrode=electrode_out;

[dummy,fstem]=fileparts(file_mat);

fprintf('\nmust load [%s] to locate electrodes in ***pre-op*** MRI!\n\n',sprintf('%s_s031.mat',fstem));
save(sprintf('%s_s031.mat',fstem),'electrode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[vol_stc,vv,d0,d1,timeVec]=inverse_read_stc(file_vol_stc);
vol_stc=etc_z(vol_stc,find(timeVec<0),'flag_baseline_correct',1);
t_idx=find((timeVec>=min(timeVec_range))&(timeVec<=max(timeVec_range)));
timeVec=timeVec(t_idx);
vol_stc=vol_stc(:,t_idx);

%sub-sampling in time
timeVec=timeVec(1:5:end);
vol_stc=vol_stc(:,1:5:end);

load(file_forward_mat);
hh={'lh','rh'};
for hemi_idx=1:2
   file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hh{hemi_idx},'orig');
   [A(hemi_idx).vertex_coords, A(hemi_idx).faces] = read_surf(file_surf);
end;

etc_render_fsbrain('overlay_stc_timeVec',timeVec,'talxfm',(talxfm),'subject',subject,'hemi',hemi,'vol',vol,'talxfm',(talxfm),'overlay_vol_stc',vol_stc,'vol_A',A,'overlay_threshold',[50 200],'overlay_smooth',5,'view_angle',[90,-20],'camposition_l',1500,'pt',[156 87 153]);

