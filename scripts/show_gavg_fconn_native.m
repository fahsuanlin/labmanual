close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri_memory/subjects');

subject='s006';

f=MRIread('hippo_fconn_native_vol_aseg_120423_gavg_hippo_left-anat.nii');
etc_render_fsbrain('subject',subject,'overlay_vol',f,'overlay_threshold',[2 5],'overlay_truncate_neg',1);

%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
view(100,15);

mni=[0 62 4 1]'; %<----target for mPFC!!
global etc_render_fsbrain;
click_vertex_vox=inv(etc_render_fsbrain.vol.vox2ras)*inv(etc_render_fsbrain.vol_pre_xfm)*inv(etc_render_fsbrain.talxfm)*mni;
click_vertex_vox=click_vertex_vox(1:3)';
        
surface_coord=etc_render_fsbrain.vol.tkrvox2ras*[click_vertex_vox(:); 1];
surface_coord=surface_coord(1:3);

vv=etc_render_fsbrain.orig_vertex_coords;
dist=sqrt(sum((vv-repmat([surface_coord(1),surface_coord(2),surface_coord(3)],[size(vv,1),1])).^2,2));
[min_dist,min_dist_idx]=min(dist);

mni_vertex=etc_render_fsbrain.orig_vertex_coords(min_dist_idx,:);
h=plot3(mni_vertex(1),mni_vertex(2),mni_vertex(3));
set(h,'markersize',44,'color','b','marker','.');
hgexport(gcf,sprintf('s006_hippo_fconn_native_vol_aseg_120423_gavg_hippo_left_mPFC.png'), hgexport('factorystyle'),'Format','png');
delete(h);
%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
view(-80,15);

mni=[-42 -64 48 1]'; %<----target for PPC!!
global etc_render_fsbrain;
click_vertex_vox=inv(etc_render_fsbrain.vol.vox2ras)*inv(etc_render_fsbrain.vol_pre_xfm)*inv(etc_render_fsbrain.talxfm)*mni;
click_vertex_vox=click_vertex_vox(1:3)';
        
surface_coord=etc_render_fsbrain.vol.tkrvox2ras*[click_vertex_vox(:); 1];
surface_coord=surface_coord(1:3);

vv=etc_render_fsbrain.orig_vertex_coords;
dist=sqrt(sum((vv-repmat([surface_coord(1),surface_coord(2),surface_coord(3)],[size(vv,1),1])).^2,2));
[min_dist,min_dist_idx]=min(dist);

mni_vertex=etc_render_fsbrain.orig_vertex_coords(min_dist_idx,:);
h=plot3(mni_vertex(1),mni_vertex(2),mni_vertex(3));
set(h,'markersize',44,'color','b','marker','.');
hgexport(gcf,sprintf('s006_hippo_fconn_native_vol_aseg_120423_gavg_hippo_left_PPC.png'), hgexport('factorystyle'),'Format','png');
delete(h);
%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
close;
clear global etc_render_fsbrain;

f=MRIread('hippo_fconn_native_vol_aseg_120423_hippo_left-anat.nii');
etc_render_fsbrain('subject',subject,'overlay_vol',f,'overlay_threshold',[2 5],'overlay_truncate_neg',1);

%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
view(100,15);

mni=[0 62 4 1]'; %<----target for mPFC!!
global etc_render_fsbrain;
click_vertex_vox=inv(etc_render_fsbrain.vol.vox2ras)*inv(etc_render_fsbrain.vol_pre_xfm)*inv(etc_render_fsbrain.talxfm)*mni;
click_vertex_vox=click_vertex_vox(1:3)';
        
surface_coord=etc_render_fsbrain.vol.tkrvox2ras*[click_vertex_vox(:); 1];
surface_coord=surface_coord(1:3);

vv=etc_render_fsbrain.orig_vertex_coords;
dist=sqrt(sum((vv-repmat([surface_coord(1),surface_coord(2),surface_coord(3)],[size(vv,1),1])).^2,2));
[min_dist,min_dist_idx]=min(dist);

mni_vertex=etc_render_fsbrain.orig_vertex_coords(min_dist_idx,:);
h=plot3(mni_vertex(1),mni_vertex(2),mni_vertex(3));
set(h,'markersize',44,'color','b','marker','.');
hgexport(gcf,sprintf('s006_hippo_fconn_native_vol_aseg_120423_hippo_left_mPFC.png'), hgexport('factorystyle'),'Format','png');
delete(h);
%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
view(-80,15);

mni=[-42 -64 48 1]'; %<----target for PPC!!
global etc_render_fsbrain;
click_vertex_vox=inv(etc_render_fsbrain.vol.vox2ras)*inv(etc_render_fsbrain.vol_pre_xfm)*inv(etc_render_fsbrain.talxfm)*mni;
click_vertex_vox=click_vertex_vox(1:3)';
        
surface_coord=etc_render_fsbrain.vol.tkrvox2ras*[click_vertex_vox(:); 1];
surface_coord=surface_coord(1:3);

vv=etc_render_fsbrain.orig_vertex_coords;
dist=sqrt(sum((vv-repmat([surface_coord(1),surface_coord(2),surface_coord(3)],[size(vv,1),1])).^2,2));
[min_dist,min_dist_idx]=min(dist);

mni_vertex=etc_render_fsbrain.orig_vertex_coords(min_dist_idx,:);
h=plot3(mni_vertex(1),mni_vertex(2),mni_vertex(3));
set(h,'markersize',44,'color','b','marker','.');
hgexport(gcf,sprintf('s006_hippo_fconn_native_vol_aseg_120423_hippo_left_PPC.png'), hgexport('factorystyle'),'Format','png');
delete(h);
%%%%%%%%%%%%%%%%%%%%% show target given MNI coordinates %%%%%%%%%%%%%%%%%%%%
close;
clear global etc_render_fsbrain;
