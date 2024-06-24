close all; clear all;

subjects_dir='/Users/fhlin/workspace/eegmri_memory/subjects';
%subjects_dir='/space_lin2/fhlin/eegmri_memory/subjects';

subject='s002';
hemi='lh';

file_surf={
    {'outer_skin.surf'},        % scalp
    {'outer_skull.surf'},       % outer skull
    {'inner_skull.surf'},       % inner skull
    {'../surf/lh.orig'},        % brain
    {'../surf/rh.orig'},        % brain
    };

output_file_surf={
    'skin',
    'skull',
    'csf',
    'gm_lh',
    'gm_rh',
    };

head_surf_stem='skin'; %file stem for scalp; must be one of output_file_surf

tissue_name={
    'Skin';     %scalp
    'Skull';    %outer skull
    'CSF';      %inner skull
    'GM_LH';    %brain
    'GM_RH';    %brain
    };

% Electric field calculations in brain stimulation based on finite elements: An optimized processing pipeline for the generation and usage of accurate individual head models Hum. Brain Mapp., 34 (4) (2013), pp. 923-935, 10.1002/hbm.21479
tissue_conductivity=[
    0.465; % scalp (below scalp)
    0.010; % skull (below outer skull)
    1.654; % CSF (below inner skull)
    0.275; % brain (below brain)
    0.275; % brain (below brain)
    ];

tissue_enclosing={
    'FreeSpace';    % ouside scalp
    'Skin';         % outside outer skull
    'Skull';        % outside inner skull
    'CSF';          % otuside brain
    'CSF';          % otuside brain
    };


file_bem='tissue_index_bem.txt';

tms_coil_name='MagVenture_MRiB91';

target_coord=[-5.1 56.7 59.4];

flag_nav=0; %open navigation window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setenv('SUBJECTS_DIR',subjects_dir);

switch(lower(hemi))
    case 'lh'
        [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/lh.orig',subjects_dir,subject));
    case 'rh'
        [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/rh.orig',subjects_dir,subject));
end;

%prepare BEM
[status, bem_obj]=etc_tms_prepare_bem(subject,file_surf,output_file_surf,tissue_name,tissue_conductivity,tissue_enclosing,file_bem);

[a,head_surf_idx]=ismember(head_surf_stem, output_file_surf);
if(a<eps)
    str=sprintf('%s,',output_file_surf{:});
    fprintf('Error! Scalp surface [%s] is not among surfaces {%s}\n', head_surf_stem, str(1:end-1));
    return;
else
    fprintf('Scalp surface [%s] found: <%03d>::<%s>\n', head_surf_stem, head_surf_idx, output_file_surf{head_surf_idx});
end;

[status, bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, file_mesh, file_meshp] =etc_tms_efield_prep_model(file_bem);

%prepare TMS coil (initialization)
[status, strcoil, coil_P, coil_t, coil_tind]=etc_tms_prepare_coil(tms_coil_name);


%%% this is preparing navigation window objects for subsequent rendering, if any
if(flag_nav)
    global etc_render_fsbrain
    if(isfield(etc_render_fsbrain,'app_tms_nav'))
        [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav(etc_render_fsbrain.app_tms_nav, strcoil, coil_P, coil_t);
    else
        [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav([], strcoil, coil_P, coil_t,'subject',subject);
    end;
else
    tms_coil_origin=[0 0 0];
    tms_coil_axis=[0 0 -1];
    tms_coil_up=[0 1 0];
    tms_coil_xfm=eye(4);
end;

% move TMS coil to target
[tms_coil_xfm_moved]=etc_tms_target_xfm_goto(target_coord, bem_obj(head_surf_idx), tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm);


%move strcoil
strcoil = etc_tms_target_xfm_apply(strcoil, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm, tms_coil_xfm_moved);
%%% this is updating navigation window objects for subsequent rendering, if any
if(flag_nav)
    results = etc_tms_target_xfm_apply_nav(etc_render_fsbrain.app_tms_nav, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2));
end;

tms_coil_center_now=tms_coil_xfm_moved(1:3,4);
tms_coil_axis_now= -tms_coil_xfm_moved(1:3,3);
tms_coil_up_now=tms_coil_xfm_moved(1:3,2);

%calculate the e-field
coords=vertex_coords./1e3;

[status, efield1]=etc_tms_efield_surf(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords,'tissue_to_plot','GM_LH');

if(flag_nav)
    %the following is for rendering E-field
    etc_render_fsbrain.overlay_vertex=efield1.vertices;
    etc_render_fsbrain.overlay_value=efield1.Etotal;

    etc_render_fsbrain.overlay_smooth=[];
    etc_render_fsbrain.overlay_source=1;
    if(isempty(etc_render_fsbrain.overlay_vol_stc))
        etc_render_fsbrain_overlay_vol_init;
    else
        etc_render_fsbrain.overlay_vol_stc=etc_render_fsbrain.overlay_value(:);
    end;

    etc_render_fsbrain.overlay_flag_render=1;
    etc_render_fsbrain_handle('update_overlay_vol');
    etc_render_fsbrain_handle('redraw');
    etc_render_fsbrain_handle('draw_pointer');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% tuning loop, if any, starts here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rotation; 30 degrees
[tms_coil_xfm_moved_tuned]=etc_tms_target_xfm_tune(target_coord, bem_obj(head_surf_idx), tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm_moved, 4, 30);

%move strcoil
strcoil = etc_tms_target_xfm_apply(strcoil, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm_moved, tms_coil_xfm_moved_tuned);
%%% this is updating navigation window objects for subsequent rendering, if any
if(flag_nav)
        results = etc_tms_target_xfm_apply_nav(etc_render_fsbrain.app_tms_nav, tms_coil_xfm_moved_tuned(1:3,4).*1e3, -tms_coil_xfm_moved_tuned(1:3,3), tms_coil_xfm_moved_tuned(1:3,2));
end;


%calculate the e-field
coords=vertex_coords./1e3;

[status, efield2]=etc_tms_efield_surf(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords,'tissue_to_plot','GM_LH');

if(flag_nav)
    %the following is for rendering E-field
    etc_render_fsbrain.overlay_vertex=efield2.vertices;
    etc_render_fsbrain.overlay_value=efield2.Etotal;

    etc_render_fsbrain.overlay_smooth=[];
    etc_render_fsbrain.overlay_source=1;
    if(isempty(etc_render_fsbrain.overlay_vol_stc))
        etc_render_fsbrain_overlay_vol_init;
    else
        etc_render_fsbrain.overlay_vol_stc=etc_render_fsbrain.overlay_value(:);
    end;

    etc_render_fsbrain.overlay_flag_render=1;
    etc_render_fsbrain_handle('update_overlay_vol');
    etc_render_fsbrain_handle('redraw');
    etc_render_fsbrain_handle('draw_pointer');
end;

save tms_script_test.mat efield1 efield2