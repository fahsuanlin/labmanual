close all; clear all;

subjects_dir='/Users/fhlin/workspace/eegmri_memory/subjects';
subjects_dir='/space_lin2/fhlin/eegmri_memory/subjects';
root_dir='/Users/fhlin/workspace/eegmri_memory/s012/tms_e_field';
root_dir='/space_lin2/fhlin/eegmri_memory/s012/tms_e_field';

subject={'s012'};
hemi={'lh','rh'};

file_surf={
    {'outer_skin.surf'},        % scalp
    {'outer_skull.surf'},       % outer skull
    {'inner_skull.surf'},       % inner skull
    {'../surf/lh.orig'},        % brain
    {'../surf/rh.orig'},        % brain
    };

%surf_category=categorical({'skull','skin','GM','WM','CSF','unspecified'});
surf_category={
    'skin';    % ouside scalp
    'skull';         % outside outer skull
    'CSF';        % outside inner skull
    'GM';          % otuside brain
    'GM';          % otuside brain
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

output_stem='mri_b91_mpfc_fc1_bh_proj_efield';

target_coord=[-13 56 35];
targeg_coord_mni=[-15 58 38]; %fc target

flag_nav=1; %open navigation window
bem_prep_status=0;
bem_def_status=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mni=targeg_coord_mni(:);
        vol=MRIread(sprintf('%s/%s/mri/orig.mgz',subjects_dir,subject{1}));
        talxfm=etc_read_xfm('subjects_dir',subjects_dir,'subject',subject{1});
        vol_pre_xfm=eye(4);
        click_vertex_vox=inv(vol.vox2ras)*inv(vol_pre_xfm)*inv(talxfm)*cat(1,mni(:),1);
        click_vertex_vox=click_vertex_vox(1:3)';
        
        surface_coord=vol.tkrvox2ras*[click_vertex_vox(:); 1];
        surface_coord=surface_coord(1:3);
        
        target_coord=surface_coord(:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setenv('SUBJECTS_DIR',subjects_dir);

for subj_idx=1:length(subject)
    if(iscell(hemi))
        v=[];
        f=[];
        for hemi_idx=1:length(hemi)
            switch(lower(hemi{hemi_idx}))
                case 'lh'
                    [vertex_tmp{hemi_idx}, face_tmp{hemi_idx}] = read_surf(sprintf('%s/%s/surf/lh.orig',subjects_dir,subject{subj_idx}));
                case 'rh'
                    [vertex_tmp{hemi_idx}, face_tmp{hemi_idx}] = read_surf(sprintf('%s/%s/surf/rh.orig',subjects_dir,subject{subj_idx}));
            end;
            f=cat(1,f,face_tmp{hemi_idx}+size(v,1));
            v=cat(1,v,vertex_tmp{hemi_idx});
        end;
        vertex_coords=v;
        faces=f;
    else
        switch(lower(hemi))
            case 'lh'
                [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/lh.orig',subjects_dir,subject{subj_idx}));
            case 'rh'
                [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/rh.orig',subjects_dir,subject{subj_idx}));
        end;
    end;

    %prepare BEM
    %path_bem=sprintf('%s/%s/analysis',root_dir,subject{subj_idx});
    path_bem=root_dir;
    [bem_def_status, bem_obj]=etc_tms_prepare_bem(subject{subj_idx},file_surf,output_file_surf,tissue_name,tissue_conductivity,tissue_enclosing,file_bem,'path_bem',path_bem,'surf_category',surf_category);

    [a,head_surf_idx]=ismember(head_surf_stem, output_file_surf);
    if(a<eps)
        str=sprintf('%s,',output_file_surf{:});
        fprintf('Error! Scalp surface [%s] is not among surfaces {%s}\n', head_surf_stem, str(1:end-1));
        return;
    else
        fprintf('Scalp surface [%s] found: <%03d>::<%s>\n', head_surf_stem, head_surf_idx, output_file_surf{head_surf_idx});
    end;


    
    pdir=pwd;
    cd(path_bem);
    [bem_prep_status, bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, file_mesh, file_meshp] =etc_tms_efield_prep_model(file_bem,'path_tissue_mesh',path_bem);
    cd(pdir);
    %prepare TMS coil (initialization)
    [status, strcoil, coil_P, coil_t, coil_tind]=etc_tms_prepare_coil(tms_coil_name);

    %%% this is preparing navigation window objects for subsequent rendering, if any
    if(flag_nav)
        global etc_render_fsbrain

        if(bem_def_status|bem_prep_status)
            %move the scalp surface to the last....
            bem_idx=[1:length(bem_obj)];
            bem_idx(head_surf_idx)=[];
            bem_idx(end+1)=head_surf_idx;
            bem_obj_tmp=bem_obj([bem_idx]);

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
        end;

        if(isfield(etc_render_fsbrain,'app_tms_nav'))
            [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav(etc_render_fsbrain.app_tms_nav, strcoil, coil_P, coil_t,'target_coord',target_coord,'hemi',hemi);
        else
            [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav([], strcoil, coil_P, coil_t,'subject',subject{subj_idx},'target_coord',target_coord,'hemi',hemi);
        end;
        if(bem_def_status|bem_prep_status)
            %etc_render_fsbrain.surf_obj=bem_obj;
            if(isfield(etc_render_fsbrain,'app_tms_nav'))
                if(~isempty(etc_render_fsbrain.app_tms_nav))
                    if(isvalid(etc_render_fsbrain.app_tms_nav))
                        if(bem_def_status)
                            etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.SurfaceDropDown));
                            etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.DefModelLamp),'g','r');
                        end;
                        if(bem_prep_status)
                            etc_render_fsbrain_tms_nav_notify(etc_render_fsbrain.app_tms_nav,struct('Source', etc_render_fsbrain.app_tms_nav.PrepModelLamp),'g','');
                        end;
                    end;
                end;
            end;
        end;
    else
        tms_coil_origin=[0 0 0];
        tms_coil_axis=[0 0 -1];
        tms_coil_up=[0 1 0];
        tms_coil_xfm=eye(4);
    end;

    % move TMS coil to target
    %[tms_coil_xfm_moved, tms_coil_xfm_mm, coil_center, coil_orientation]=etc_tms_target_xfm_goto(target_coord, bem_obj(head_surf_idx), tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm,'coil_center',surf_target_coord);
    [tms_coil_xfm_moved, tms_coil_xfm_mm, coil_center, coil_orientation]=etc_tms_target_xfm_goto(target_coord, bem_obj(head_surf_idx), tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm);


    %move strcoil
    [strcoil, strcoil_xfm] = etc_tms_target_xfm_apply(strcoil, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm, tms_coil_xfm_moved);
    %%% this is updating navigation window objects for subsequent rendering, if any
    if(flag_nav)
        results = etc_tms_target_xfm_apply_nav(etc_render_fsbrain.app_tms_nav, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2));
    end;

    tms_coil_center_now=tms_coil_xfm_moved(1:3,4);
    tms_coil_axis_now= -tms_coil_xfm_moved(1:3,3);
    tms_coil_up_now=tms_coil_xfm_moved(1:3,2);

    %calculate the e-field
    coords=vertex_coords./1e3;

    %[efield_status, efield1]=etc_tms_efield_surf(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords,'tissue_to_plot','GM_LH');
    [efield_status, efield1]=etc_tms_efield_surf(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords,'tissue_to_plot',{'GM_LH','GM_RH'});

    if(iscell(hemi))
%         switch(lower(hemi{hemi_idx}))
%             case 'lh'
%                 efield.lh=efield1;
%             case 'rh'
%                 efield.rh=efield1;
%         end;
        efield=efield1;
    else
        efield=efield1;
    end;

    %append TMS coil object
    efield.tms_coil_xfm=tms_coil_xfm_moved;
    efield.tms_coil_name=tms_coil_name;
    
    %append TMS target coord.
    efield.target_coord=target_coord;

    %append BEM
    efield.bem_obj=bem_obj;
    efield.scalp_bem_index=head_surf_idx;
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

        axis(etc_render_fsbrain.brain_axis,'tight');
        xlim=get(etc_render_fsbrain.brain_axis,'xlim');
        ylim=get(etc_render_fsbrain.brain_axis,'ylim');
        zlim=get(etc_render_fsbrain.brain_axis,'zlim');

        hgexport(gcf,sprintf('%s_init.png',output_stem), hgexport('factorystyle'),'Format','png');
    end;

end;

pdir=pwd;
fn=sprintf('%s/%s.mat',root_dir,output_stem);
fprintf('saving [%s]...\n',fn);
save(fn, 'efield');
cd(pdir);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rotation_deg=[0:15:360-1];
     
for rot_idx=1:length(rotation_deg)
    if(flag_nav)
        [tms_coil_xfm_moved_tuned]=etc_tms_target_xfm_tune(target_coord, bem_obj(head_surf_idx), tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm_moved, 4, rotation_deg(rot_idx),'app',etc_render_fsbrain.app_tms_nav);
    else
        [tms_coil_xfm_moved_tuned]=etc_tms_target_xfm_tune(target_coord, bem_obj(head_surf_idx), tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm_moved, 4, rotation_deg(rot_idx));
    end;

    %calculate the e-field
    coords=vertex_coords./1e3;

    [status, efield1]=etc_tms_efield_surf(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords,'tissue_to_plot',{'GM_LH','GM_RH'});
    
    efield_rot{rot_idx}=efield1;

    %append TMS coil object
    efield_rot{rot_idx}.tms_coil_xfm=tms_coil_xfm_moved_tuned;
    efield_rot{rot_idx}.tms_coil_name=tms_coil_name;
    
    %append TMS target coord.
    efield_rot{rot_idx}.target_coord=target_coord;

    %append BEM
    efield_rot{rot_idx}.bem_obj=bem_obj;
    efield_rot{rot_idx}.scalp_bem_index=head_surf_idx;


    if(flag_nav)
        %the following is for rendering E-field
        etc_render_fsbrain.overlay_vertex=efield_rot{rot_idx}.vertices;
        etc_render_fsbrain.overlay_value=efield_rot{rot_idx}.Etotal;

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

                set(etc_render_fsbrain.brain_axis,'xlim',xlim);
        set(etc_render_fsbrain.brain_axis,'ylim',ylim);
        set(etc_render_fsbrain.brain_axis,'zlim',zlim);

        rot_str=sprintf('rot%03d',rotation_deg(rot_idx));

        hgexport(gcf,sprintf('%s_%s.png',output_stem, rot_str), hgexport('factorystyle'),'Format','png');

    end;
end;

pdir=pwd;
fn=sprintf('%s/%s.mat',root_dir,output_stem);
fprintf('saving [%s]...\n',fn);
save(fn,'-append','efield_rot','rotation_deg');
cd(pdir);
