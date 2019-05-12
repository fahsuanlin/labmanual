clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/sinica_meg/subjects/'); %for MAC/Linux

%source space
file_source_fif='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/041019-5-src.fif';

surfin_lh='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/surf/lh.white';
surfin_rh='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/surf/rh.white';
surf_outer_skin='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/watershed/041019_outer_skin_surface';
surf_outer_skull='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/watershed/041019_outer_skull_surface';
surf_inner_skull='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/watershed/041019_inner_skull_surface';


eeg_output_name='eeg_fwd_prep_042319.mat';


h_mri_brain_lh=[];
h_mri_brain_rh=[];
h_mri_head=[];
h_mri_source_location=[];
h_mri_source_norm=[];
h_meg_array=[];


flag_show_outer_skin=0;
flag_show_outer_skull=1;
flag_show_inner_skull=1;

flag_show_mri_brain=1;
flag_show_mri_head=1;
flag_show_meg_array=0;
flag_show_meg_source_location=1;
flag_show_meg_source_norm=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[src] = mne_read_source_spaces(file_source_fif);

%desired_ntri=2000;
[verts_lh, faces_lh] = mne_read_surface(surfin_lh);
[verts_rh, faces_rh] = mne_read_surface(surfin_rh);
[verts_osc, faces_osc] = mne_read_surface(surf_outer_skin);
[verts_osk, faces_osk] = mne_read_surface(surf_outer_skull);
[verts_isk, faces_isk] = mne_read_surface(surf_inner_skull);



if(flag_show_outer_skin)
    h_outer_skin=patch('Faces',faces_osc,'Vertices',verts_osc,'FaceColor',[1 1 1].*0.8,'facealpha',0.2);
    set(h_outer_skin,'facecolor',[256 180 100]./256,'edgecolor','none','facealpha',0.2);
end;

if(flag_show_outer_skull)
    h_outer_skull=patch('Faces',faces_osk,'Vertices',verts_osk,'FaceColor',[1 1 1].*0.8,'facealpha',0.2);
    set(h_outer_skull,'facecolor',[256 180 100]./256,'edgecolor','none','facealpha',0.2);
end;

if(flag_show_inner_skull)
    h_inner_skull=patch('Faces',faces_isk,'Vertices',verts_isk,'FaceColor',[1 1 1].*0.8,'facealpha',0.2);
    set(h_inner_skull,'facecolor',[256 180 100]./256,'edgecolor','none','facealpha',1);
end;
axis equal off vis3d; camlight; material dull;

load(eeg_output_name);
elec.elecpos=points;
elec.label=points_label;


etc_render_topo('vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_aux_point_coords',elec.elecpos,'topo_aux_point_name',elec.label);


return;
