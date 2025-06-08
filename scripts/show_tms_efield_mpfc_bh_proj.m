close all; clear all;

subject='s003';
hemi={'lh','rh'};

tms_coil_name='MagVenture_MRiB91';


%%%%%%%%%%%%%%%%%%%
load mri_b91_mpfc_bh_proj_efield.mat; %variable 'efield' must be there


%%%%%%%%%%%%%%%%%%%
% BEM
global etc_render_fsbrain

bem_idx=[1:length(efield.bem_obj)];
bem_idx(efield.scalp_bem_index)=[];
bem_idx(end+1)=efield.scalp_bem_index;
bem_obj_tmp=efield.bem_obj([bem_idx]);

for idx=1:length(bem_obj_tmp)
    [pp,ff]=fileparts(bem_obj_tmp(idx).filename);
    app.SurfaceDropDown.Items{idx}=ff;
    app.SurfaceDropDown.Value=ff;

    [vv,ff]=read_surf(bem_obj_tmp(idx).filename);

    cc=colororder;
    cc_idx=length(bem_obj_tmp);
    cc_idx=mod(cc_idx+1,size(cc,1));
    if(cc_idx==0) cc_idx=size(cc,1); end;
    bem_obj_tmp(idx).color=cc(cc_idx,:);

    if(idx==length(bem_obj_tmp))
        bem_obj_tmp(idx).flag_show=1;
    else
        bem_obj_tmp(idx).flag_show=0;
    end;

    bem_obj_tmp(idx).patch_handle=patch('vertices',vv,'faces',ff+1,'edgecolor','none','facecolor',bem_obj_tmp(idx).color,'facealpha',0.2,'visible',bem_obj_tmp(idx).flag_show); %show scalp
end;

etc_render_fsbrain.surf_obj=bem_obj_tmp;

%etc_render_fsbrain_init('subject',subject,'hemi',hemi);
%%%%%%%%%%%%%%%%%%%
% tms coil object

[status, strcoil, coil_P, coil_t, coil_tind]=etc_tms_prepare_coil(efield.tms_coil_name,'flag_save',0);


%[status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm, strcoil] = etc_tms_show_tms_coil(coil_P, coil_t,'target_coord',efield.target_coord,'hemi',hemi,'tms_coil_xfm',efield.tms_coil_xfm,'strcoil',strcoil);
[status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav([], strcoil, coil_P, coil_t,'subject',subject,'target_coord',efield.target_coord,'hemi',hemi);

global etc_render_fsbrain;
if(isfield(etc_render_fsbrain,'app_tms_nav'))
    if(~isempty(etc_render_fsbrain.app_tms_nav))
        if(isvalid(etc_render_fsbrain.app_tms_nav))etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.SurfaceDropDown));
            etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.DefModelLamp),'g','r');
            etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.PrepModelLamp),'g','');
        end;
    end;
end;

results = etc_tms_target_xfm_apply_nav(etc_render_fsbrain.app_tms_nav, efield.tms_coil_xfm(1:3,4).*1e3, -efield.tms_coil_xfm(1:3,3), efield.tms_coil_xfm(1:3,2));


%%%%%%%%%%%%%%%%%%%
% show e-field

global etc_render_fsbrain;
%the following is for rendering E-field
etc_render_fsbrain.overlay_vertex=efield.vertices;
etc_render_fsbrain.overlay_value=efield.Etotal;

etc_render_fsbrain.overlay_smooth=[];
etc_render_fsbrain.overlay_source=1;

etc_render_fsbrain_overlay_vol_init;

etc_render_fsbrain.overlay_flag_render=1;
etc_render_fsbrain_handle('update_overlay_vol');
etc_render_fsbrain_handle('redraw');
etc_render_fsbrain_handle('draw_pointer');

axis(etc_render_fsbrain.brain_axis,'tight');
xlim=get(etc_render_fsbrain.brain_axis,'xlim');
ylim=get(etc_render_fsbrain.brain_axis,'ylim');
zlim=get(etc_render_fsbrain.brain_axis,'zlim');

%%%%%%%%%%%%%%%%%%%

