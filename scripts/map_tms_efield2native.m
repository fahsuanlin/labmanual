close all; clear all;

subject='s012';
hemi={'lh','rh'};
file_efield={
    'mri_b91_mpfc_atlas_bh_proj_efield.mat';
    'mri_b91_mpfc_fc1_bh_proj_efield.mat';
    'mri_b91_mpfc_fc2_bh_proj_efield.mat';
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for file_idx=1:length(file_efield)
    load(file_efield{file_idx});

    brain_fig=figure; set(brain_fig,'visible','off');
    etc_render_fsbrain('subject',subject,'hemi',hemi);


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

    [dummy,fstem]=fileparts(file_efield{file_idx});
    fn=sprintf(sprintf('%s_anat.nii',fstem))
    MRIwrite(etc_render_fsbrain.overlay_vol,fn);

    if(exist('efield_rot','var'))
        vol_rot=etc_render_fsbrain.overlay_vol;
        for rot_idx=1:length(efield_rot)
            etc_render_fsbrain.overlay_vertex=efield_rot{rot_idx}.vertices;
            etc_render_fsbrain.overlay_value=efield_rot{rot_idx}.Etotal;

            %etc_render_fsbrain_overlay_vol_init;
            etc_render_fsbrain.overlay_vol_stc=efield_rot{rot_idx}.Etotal;

            etc_render_fsbrain.overlay_flag_render=1;
            etc_render_fsbrain_handle('update_overlay_vol');
            etc_render_fsbrain_handle('redraw');
            etc_render_fsbrain_handle('draw_pointer');
            vol_rot.vol=cat(4,vol_rot.vol,etc_render_fsbrain.overlay_vol.vol);
            vol_rot.nframes=size(vol_rot.vol,4);
        end;
        fn=sprintf(sprintf('%s_rot_anat.nii',fstem));
        MRIwrite(vol_rot,fn);

    end;
    clear global etc_render_fsbrain;
end;

%%%%%%%%%%%%%%%%%%%

