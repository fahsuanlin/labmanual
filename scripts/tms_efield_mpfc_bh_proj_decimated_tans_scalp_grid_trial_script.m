close all; clear all;

root_dir=fileparts(mfilename('fullpath'));
if(isempty(root_dir))
    root_dir=pwd;
end;
subjects_dir=getenv('SUBJECTS_DIR');
if(isempty(subjects_dir))
    subjects_dir=fullfile(fileparts(fileparts(root_dir)),'subjects');
end;

subject_dir=fileparts(root_dir);
[~,subject_name]=fileparts(subject_dir);
subject={subject_name};
hemi={'lh','rh'};

file_surf={
    {'outer_skin.surf'},        % scalp
    {'outer_skull.surf'},       % outer skull
    {'inner_skull.surf'},       % inner skull
    {'../surf/lh.orig'},        % brain
    {'../surf/rh.orig'},        % brain
    };

surf_category={
    'skin';          % outside scalp
    'skull';         % outside outer skull
    'CSF';           % outside inner skull
    'GM';            % outside brain
    'GM';            % outside brain
    };

output_file_surf={
    'skin',
    'skull',
    'csf',
    'gm_lh',
    'gm_rh',
    };

head_surf_stem='skin'; % file stem for scalp; must be one of output_file_surf

tissue_name={
    'Skin';     % scalp
    'Skull';    % outer skull
    'CSF';      % inner skull
    'GM_LH';    % brain
    'GM_RH';    % brain
    };

% Electric field calculations in brain stimulation based on finite elements:
% An optimized processing pipeline for the generation and usage of accurate
% individual head models. Hum. Brain Mapp., 34 (4) (2013), pp. 923-935.
tissue_conductivity=[
    0.465; % scalp (below scalp)
    0.010; % skull (below outer skull)
    1.654; % CSF (below inner skull)
    0.275; % brain (below brain)
    0.275; % brain (below brain)
    ];

tissue_enclosing={
    'FreeSpace';    % outside scalp
    'Skin';         % outside outer skull
    'Skull';        % outside inner skull
    'CSF';          % outside brain
    'CSF';          % outside brain
    };

file_bem='tissue_index_bem.txt';
tms_coil_name='MagVenture_MRiB91';

output_stem='tms_efield_mpfc_bh_proj_decimated_dfconn_sup_mPFC_tans_scalp_grid_trial';
flag_cleanup_previous_optimize_outputs=1;

% Functional ROI definition:
%   compare_dfconn_hippo_entorhinal_060525-anat.nii > volume_threshold
%   intersected with average/sup_mPFC_hippo_ent-?h.label morphed to s025.
average_dir=fullfile(fileparts(fileparts(root_dir)),'average');
batch_log_dir=fullfile(average_dir,'batch_logs');
if(~exist(batch_log_dir,'dir'))
    mkdir(batch_log_dir);
end;
fsaverage_label_lh=fullfile(average_dir,'sup_mPFC_hippo_ent-lh.label');
fsaverage_label_rh=fullfile(average_dir,'sup_mPFC_hippo_ent-rh.label');
dfconn_volume_file=fullfile(fileparts(root_dir),'resting_analysis','compare_dfconn_hippo_entorhinal_060525-anat.nii');
volume_threshold=1.0;
volume_surface_name='orig';
volume_projfrac=0;
volume_interp='nearest';
target_label_prefix='compare_dfconn_hippo_entorhinal_060525_gt1_sup_mPFC_hippo_ent';

% Atlas baseline comparison. This follows the existing atlas script convention:
% MNI305 target coordinates are converted to the subject native surface
% coordinate before the coil is placed.
flag_include_atlas_comparison=0;
atlas_target_coord_mni=[0 62 21];
atlas_coil_rotation_deg=0;

% TANS scalp search grid: project the functional-ROI centroid to the scalp,
% then sample a disk on outer_skin.surf.  These are scalp locations, not
% cortical ROI vertices projected independently to the scalp.
scalp_search_radius_mm=20;
scalp_grid_spacing_mm=2;
candidate_coil_rotation_deg=-90:15:90;

% Two orthogonal basis poses require BEM solves; all requested roll angles
% are synthesized from their vector E-fields.
flag_two_basis_orientation_synthesis=1;
flag_write_search_checkpoint=1;
search_checkpoint_interval=1; % save after this many unique ROI positions

% Exact-vs-synthesized validation on the first few scalp positions.
flag_validate_two_basis_synthesis=1;
validation_candidate_count=10;

top_result_count=5;
flag_write_candidate_list=1;
target_coord=[];
target_coord_mni=[];
target_coil_rotation_deg=[];
target_candidate_info=[];

% These fractions reduce the BEM surfaces before the charge solve. This is
% the part that controls runtime in etc_tms_efield_surf.
head_surface_reduce_fraction=0.25;  % scalp/skull/CSF surfaces
brain_surface_reduce_fraction=0.10; % lh/rh GM surfaces

% Reuse the reduced CombinedMesh files if they already exist for this
% reduction level. Set these to 1 after changing surfaces or reduce fractions.
flag_force_rebuild_decimated_surfaces=0;
flag_force_rebuild_prepared_model=0;

% If true, efield is expanded back to all original brain vertices by nearest
% sparse sample followed by inverse_smooth on the full surface.
flag_interpolate_to_full_brain=1;
nearest_chunk_size=2000; % fallback chunk size when knnsearch is unavailable
full_interpolation_smooth_step=5;
flag_stc=0; % keep STC export disabled during the search pass
% Native NIfTI export is an optional, expensive post-processing step: each
% retained result launches mri_surf2vol once per E-field scalar plus once for
% the ROI. Keep it off during optimization; set to 1 when these volumes are
% explicitly needed after the ranked MAT/CSV outputs have been produced.
flag_native_nii=0;
% If true, write four native-space target-score volumes after all candidates:
% all rotations (4D), maximum TANS, mean TANS, and best rotation (degrees).
flag_native_tans_nii=1;

flag_nav=0; % open navigation window and render the sparse e-field
overlay_threshold='tans'; % 'tans' uses [p99 p99.9] Etotal; [] uses automatic threshold

% On-target ROI calculation. If label files are provided, they define the
% ROI. Otherwise the ROI is a cortical disk centered on target_coord_mni
% after conversion to this subject's surface coordinate.
flag_on_target=1;
on_target_label_files={}; % set after dfconn>threshold superior-mPFC ROI labels are created
on_target_circle_diameter_mm=30; % 30 mm diameter
on_target_circle_hemi='both'; % 'both', 'auto', 'lh', or 'rh'
on_target_exclude_sulcal_vertices=0;
on_target_sulc_threshold=0;
on_target_sulc_keep_direction='negative'; % FreeSurfer .sulc: keep <=0 for gyral crown; use 'positive' for TANS CIFTI convention
on_target_empty_after_sulcal_mask_action='nearest_gyral'; % 'nearest_gyral' or 'error' for circle ROIs

% TANS on-target score: surface-area weighted overlap between the E-field
% magnitude hotspot and the target ROI, averaged over 99.9%..99.0% cutoffs.
on_target_field_name='Etotal'; % TANS uses E-field magnitude (SimNIBS magnE)
on_target_hotspot_percentile_thresholds=linspace(99.9,99,10);
on_target_hotspot_hemi='both'; % 'both', 'auto', 'lh', or 'rh'
flag_on_target_write_hotspot_label=0; % search pass: do not write hotspot labels for every candidate
on_target_hotspot_label_percentile=99.5; % paper visualization threshold

% dI/dt dose curve: for a fixed coil pose, E-field scales linearly with
% dI/dt. etc_tms_efield_surf currently computes fields at 9.4e7 A/s
% (=94 A/us), so these settings estimate the absolute-threshold hotspot
% on-target percentage at each requested dI/dt without recomputing BEM.
flag_on_target_didt=0;
on_target_didt_reference_A_per_us=94;
on_target_didt_A_per_us=20:5:120;
on_target_didt_activation_threshold_V_per_m=[50 100];
on_target_didt_hotspot_hemi='both'; % 'both', 'auto', 'lh', or 'rh'

bem_prep_status=0;
bem_def_status=0;
on_target=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setenv('SUBJECTS_DIR',subjects_dir);

surface_reduce_fraction=[
    head_surface_reduce_fraction;
    head_surface_reduce_fraction;
    head_surface_reduce_fraction;
    brain_surface_reduce_fraction;
    brain_surface_reduce_fraction;
    ];

path_bem=sprintf('%s/bem_decimated_head%03d_brain%03d',root_dir,round(head_surface_reduce_fraction.*1000),round(brain_surface_reduce_fraction.*1000));

if(flag_cleanup_previous_optimize_outputs)
    etc_local_cleanup_previous_optimize_outputs(root_dir, output_stem);
end;

for subj_idx=1:length(subject)
    [bem_def_status, bem_obj, brain_decimation]=etc_local_prepare_or_load_decimated_bem( ...
        subjects_dir, subject{subj_idx}, file_surf, output_file_surf, tissue_name, ...
        tissue_conductivity, tissue_enclosing, file_bem, path_bem, surf_category, ...
        surface_reduce_fraction, flag_force_rebuild_decimated_surfaces, nearest_chunk_size);

    fprintf('Decimated BEM path: [%s]\n',path_bem);
    fprintf('Full brain mesh: %d vertices, %d faces\n', brain_decimation.full_nvertices, brain_decimation.full_nfaces);
    fprintf('Sparse brain mesh: %d vertices, %d faces\n', brain_decimation.sparse_nvertices, brain_decimation.sparse_nfaces);
    fprintf('Total BEM faces: %d full -> %d decimated\n', brain_decimation.total_full_bem_faces, brain_decimation.total_sparse_bem_faces);

    label_config=struct;
    label_config.subjects_dir=subjects_dir;
    label_config.fsaverage_lh_label=fsaverage_label_lh;
    label_config.fsaverage_rh_label=fsaverage_label_rh;
    label_config.volume_threshold=volume_threshold;
    label_config.surface_name=volume_surface_name;
    label_config.projfrac=volume_projfrac;
    label_config.vol2surf_interp=volume_interp;
    label_config.label_prefix=target_label_prefix;

    if(exist(dfconn_volume_file,'file')~=2)
        error('Cannot find dfconn volume [%s].',dfconn_volume_file);
    end;
    label_info=etc_local_create_dfconn_sup_mPFC_labels(subject{subj_idx}, root_dir, dfconn_volume_file, brain_decimation, label_config, batch_log_dir);
    if(label_info.total_nvertices<1)
        error('The dfconn>%.3g superior-mPFC ROI is empty for %s.',volume_threshold,subject{subj_idx});
    end;
    on_target_label_files=label_info.label_files(1:2);

    [a,head_surf_idx]=ismember(head_surf_stem, output_file_surf);
    if(a<eps)
        str=sprintf('%s,',output_file_surf{:});
        error('Scalp surface [%s] is not among surfaces {%s}.',head_surf_stem,str(1:end-1));
    end;
    [scalp_vertex_full,scalp_face_full]=etc_local_read_source_surface( ...
        subjects_dir, subject{subj_idx}, file_surf{head_surf_idx});

    [target_coords, target_coords_mni, target_coil_rotation_deg, target_candidate_info, native_to_mni_xfm]=etc_local_build_tans_scalp_candidates( ...
        subjects_dir, subject{subj_idx}, label_info, brain_decimation, ...
        scalp_vertex_full, scalp_face_full, scalp_search_radius_mm, ...
        scalp_grid_spacing_mm, candidate_coil_rotation_deg);
    if(flag_include_atlas_comparison)
        [target_coords, target_coords_mni, target_coil_rotation_deg, target_candidate_info]=etc_local_prepend_atlas_comparison_target( ...
            subjects_dir, subject{subj_idx}, atlas_target_coord_mni, atlas_coil_rotation_deg, ...
            target_coords, target_coords_mni, target_coil_rotation_deg, target_candidate_info);
    end;
    n_targets=size(target_coords,1);
    fprintf('Functional ROI labels: lh=%d/%d, rh=%d/%d, total=%d; dfconn mean=%1.3f, max=%1.3f.\n', ...
        label_info.hemi(1).nvertices, label_info.hemi(1).base_label_nvertices, ...
        label_info.hemi(2).nvertices, label_info.hemi(2).base_label_nvertices, ...
        label_info.total_nvertices, label_info.dfconn_summary.mean, label_info.dfconn_summary.max);
    fprintf('Evaluating %d target poses: %d atlas comparison + %d optimized poses (%d unique ROI coordinates x %d coil rotations).\n', ...
        n_targets, target_candidate_info.n_atlas_targets, target_candidate_info.n_optimized_targets, ...
        target_candidate_info.n_unique_coords, target_candidate_info.n_rotations);
    if(flag_write_candidate_list)
        etc_local_write_target_candidate_csv(fullfile(root_dir,sprintf('%s_candidates.csv',output_stem)), target_candidate_info);
    end;

    fprintf('Scalp surface [%s] found: <%03d>::<%s>; full vertices=%d.\n', ...
        head_surf_stem, head_surf_idx, output_file_surf{head_surf_idx},size(scalp_vertex_full,1));

    flag_force_rebuild_prepared_model_now=flag_force_rebuild_prepared_model || brain_decimation.flag_decimated_surfaces_rebuilt;
    [bem_prep_status, bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC]=etc_local_load_or_prepare_model(file_bem,path_bem,flag_force_rebuild_prepared_model_now);

    % Prepare TMS coil and optionally initialize the navigation UI.
    [status, strcoil, coil_P, coil_t, coil_tind]=etc_tms_prepare_coil(tms_coil_name);
    strcoil_base=strcoil;
    target_coord_nav_init=target_coords(1,:);

    if(flag_nav)
        global etc_render_fsbrain

        etc_local_init_nav_bem_surfaces(bem_obj, head_surf_idx);

        if(isfield(etc_render_fsbrain,'app_tms_nav'))
            [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav(etc_render_fsbrain.app_tms_nav, strcoil, coil_P, coil_t,'target_coord',target_coord_nav_init,'hemi',hemi);
        else
            [status, tms_coil_origin, tms_coil_axis, tms_coil_up, tms_coil_xfm] = etc_tms_init_nav([], strcoil, coil_P, coil_t,'subject',subject{subj_idx},'target_coord',target_coord_nav_init,'hemi',hemi);
        end;

        etc_local_notify_nav_model_status(bem_def_status, bem_prep_status);


        % for tms_nav window....
        %update bem info.
        etc_render_fsbrain.tissue_def_file=file_bem;
        for tissue_idx=1:length(tissue_name)
            etc_render_fsbrain.tissue_obj(tissue_idx).tissue_name=tissue_name{tissue_idx};
            etc_render_fsbrain.tissue_obj(tissue_idx).tissue_mat=bem_obj(tissue_idx).filemat;
            etc_render_fsbrain.tissue_obj(tissue_idx).tissue_cond=tissue_conductivity(tissue_idx);
            etc_render_fsbrain.tissue_obj(tissue_idx).tissue_enclosing_tissue=tissue_enclosing{tissue_idx};
        end;

    else
        tms_coil_origin=[0 0 0];
        tms_coil_axis=[0 0 -1];
        tms_coil_up=[0 1 0];
        tms_coil_xfm=eye(4);
    end;

    % Calculate the e-field on the reduced brain vertices. The expensive
    % BEM solve now uses the reduced surfaces above.
    coords=brain_decimation.sparse_vertex_coords./1e3;
    tissue_to_plot=etc_local_hemi_to_tissue(hemi);
    brain_decimation_base=brain_decimation;
    if(flag_nav || flag_interpolate_to_full_brain || flag_on_target || flag_stc)
        brain_decimation_base=etc_local_ensure_full_to_sparse_map(brain_decimation_base, nearest_chunk_size, full_interpolation_smooth_step);
    end;
    render_view_state=[];
    target_summary=repmat(etc_local_empty_target_summary_row(),0,1);
    top_result_cache={};
    basis_cache=cell(1,max(1,target_candidate_info.n_unique_coords));
    basis_position_count=0;
    basis_validation_rows=repmat(etc_local_empty_basis_validation_row(),0,1);
    search_checkpoint_file=fullfile(root_dir,sprintf('%s_two_basis_checkpoint.mat',output_stem));
    target_start_idx=1;
    checkpoint_signature=sprintf('%d_%d_%s',n_targets,target_candidate_info.n_unique_coords,mat2str(candidate_coil_rotation_deg));
    if(flag_write_search_checkpoint && exist(search_checkpoint_file,'file')==2)
        checkpoint=load(search_checkpoint_file);
        % Also accept n_targets+1: the search may have completed and then
        % failed during the final export stage.  In that case rerunning the
        % script should reuse the completed search instead of recalculating
        % every expensive E-field candidate.
        if(isfield(checkpoint,'checkpoint_signature') && strcmp(checkpoint.checkpoint_signature,checkpoint_signature) && ...
                isfield(checkpoint,'next_target_idx') && checkpoint.next_target_idx<=n_targets+1)
            target_summary=checkpoint.target_summary;
            top_result_cache=checkpoint.top_result_cache;
            basis_cache=checkpoint.basis_cache;
            target_start_idx=checkpoint.next_target_idx;
            fprintf('Resuming two-basis search from target %d/%d using checkpoint [%s].\n', ...
                target_start_idx,n_targets,search_checkpoint_file);
        end;
    end;

    for target_idx=target_start_idx:n_targets
        target_coord_this=target_coords(target_idx,:);
        target_coil_rotation_deg_this=target_coil_rotation_deg(target_idx);
        target_source_this=char(target_candidate_info.table.target_source(target_idx));
        output_stem_target=etc_local_target_output_stem(output_stem, target_idx, n_targets);
        brain_decimation=brain_decimation_base;
        on_target=[];
        efield_ontarget=[];
        brain_decimation_ontarget=[];
        stc_info=[];

        fprintf('\nTarget %d/%d: [%1.2f %1.2f %1.2f], coil rotation %1.2f deg -> [%s]\n', ...
            target_idx, n_targets, target_coord_this(1), target_coord_this(2), target_coord_this(3), target_coil_rotation_deg_this, output_stem_target);

        scalp_center_this=[ ...
            double(target_candidate_info.table.scalp_center_x(target_idx)), ...
            double(target_candidate_info.table.scalp_center_y(target_idx)), ...
            double(target_candidate_info.table.scalp_center_z(target_idx))];
        scalp_normal_this=[ ...
            double(target_candidate_info.table.scalp_normal_x(target_idx)), ...
            double(target_candidate_info.table.scalp_normal_y(target_idx)), ...
            double(target_candidate_info.table.scalp_normal_z(target_idx))];
        scalp_up_this=[ ...
            double(target_candidate_info.table.scalp_up_x(target_idx)), ...
            double(target_candidate_info.table.scalp_up_y(target_idx)), ...
            double(target_candidate_info.table.scalp_up_z(target_idx))];
        [tms_coil_xfm_moved, tms_coil_xfm_mm, coil_center, coil_orientation]= ...
            etc_local_scalp_coil_xfm_goto(scalp_center_this,scalp_normal_this, ...
            scalp_up_this,tms_coil_xfm);
        tms_coil_xfm_goto=tms_coil_xfm_moved;
        if(abs(target_coil_rotation_deg_this)>eps)
            tms_coil_xfm_moved=etc_local_tms_rotate_about_axis(tms_coil_xfm_goto, target_coil_rotation_deg_this);
            tms_coil_xfm_mm=etc_local_tms_xfm_m_to_mm(tms_coil_xfm_moved);
            coil_center=tms_coil_xfm_mm(1:3,4);
            coil_orientation=(-tms_coil_xfm_moved(1:3,3)).';
        end;

        % coil_center is the scalp-contact/origin position in subject-native
        % surface RAS (mm). Transform both coil position and orientation to
        % MNI305 for navigation/reporting alongside the target coordinates.
        coil_center_mni=etc_local_apply_native_to_mni(native_to_mni_xfm, coil_center);
        coil_orientation_mni=etc_local_apply_native_to_mni_vector(native_to_mni_xfm, coil_orientation);
        coil_up=tms_coil_xfm_moved(1:3,2).';
        coil_up_mni=etc_local_apply_native_to_mni_vector(native_to_mni_xfm, coil_up);
        brainsight_pose=etc_local_compute_brainsight_pose( ...
            target_coord_this, coil_center, coil_up);

        strcoil_target=strcoil_base;
        [strcoil_target, strcoil_xfm] = etc_tms_target_xfm_apply(strcoil_target, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2), tms_coil_xfm, tms_coil_xfm_moved);
        if(flag_nav)
            etc_local_apply_coil_to_nav(tms_coil_xfm_moved);
            etc_local_update_nav_target_coordinate(target_coord_this);
            etc_local_set_nav_efield_status('r');
            etc_local_set_nav_efield_status('y');
        end;

        tic;
        basis_solve_elapsed_sec=nan;
        candidate_id_this=double(target_candidate_info.table.candidate_id(target_idx));
        if(flag_two_basis_orientation_synthesis && candidate_id_this>0)
            if(candidate_id_this>length(basis_cache))
                basis_cache{candidate_id_this}=[];
            end;
            if(isempty(basis_cache{candidate_id_this}))
                % Basis 0: the target-derived pose before roll rotation.
                [strcoil_basis0,~]=etc_tms_target_xfm_apply( ...
                    strcoil_base, tms_coil_xfm_goto(1:3,4).*1e3, ...
                    -tms_coil_xfm_goto(1:3,3), tms_coil_xfm_goto(1:3,2), ...
                    tms_coil_xfm, tms_coil_xfm_goto);
                tms_coil_xfm_basis90=etc_local_tms_rotate_about_axis(tms_coil_xfm_goto,90);
                [strcoil_basis90,~]=etc_tms_target_xfm_apply( ...
                    strcoil_base, tms_coil_xfm_basis90(1:3,4).*1e3, ...
                    -tms_coil_xfm_basis90(1:3,3), tms_coil_xfm_basis90(1:3,2), ...
                    tms_coil_xfm, tms_coil_xfm_basis90);

                basis_tic=tic;
                [status_basis0,efield_basis0]=etc_local_run_efield_solve( ...
                    bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, ...
                    cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, ...
                    RnumberE, ineighborE, EC, coords, tissue_to_plot, strcoil_basis0);
                [status_basis90,efield_basis90]=etc_local_run_efield_solve( ...
                    bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, ...
                    cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, ...
                    RnumberE, ineighborE, EC, coords, tissue_to_plot, strcoil_basis90);
                if(~isfield(efield_basis0,'E') || ~isfield(efield_basis90,'E'))
                    error('Two-basis synthesis requires vector E-fields in efield_sparse.E.');
                end;
                basis_cache{candidate_id_this}=struct( ...
                    'efield0',efield_basis0,'efield90',efield_basis90, ...
                    'tms_coil_xfm0',tms_coil_xfm_goto, ...
                    'tms_coil_xfm90',tms_coil_xfm_basis90, ...
                    'solve_elapsed_sec',toc(basis_tic), ...
                    'status0',status_basis0,'status90',status_basis90);
                basis_position_count=basis_position_count+1;
                brain_decimation.efield_elapsed_sec=nan;
                brain_decimation.efield_basis_elapsed_sec=nan;
            end;
            basis_entry=basis_cache{candidate_id_this};
            efield_sparse=etc_local_synthesize_efield_from_basis( ...
                basis_entry.efield0,basis_entry.efield90,target_coil_rotation_deg_this);
            efield_status=1;
            brain_decimation.efield_elapsed_sec=toc;
            brain_decimation.efield_basis_elapsed_sec=basis_entry.solve_elapsed_sec;
            basis_solve_elapsed_sec=basis_entry.solve_elapsed_sec;
            brain_decimation.efield_status=efield_status;
        else
            try
                [efield_status, efield_sparse]=etc_local_run_efield_solve( ...
                    bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, ...
                    cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, ...
                    RnumberE, ineighborE, EC, coords, tissue_to_plot, strcoil_target);
            catch ME
                if(flag_nav)
                    etc_local_set_nav_efield_status('r');
                end;
                rethrow(ME);
            end;
            brain_decimation.efield_elapsed_sec=toc;
            brain_decimation.efield_status=efield_status;
        end;
        % etc_tms_efield_surf marks the lamp green internally on return, but
        % the navigation E-field fields are not installed until below.
        if(flag_nav)
            etc_local_set_nav_efield_status('y');
        end;
        fprintf('Sparse E-field calculation took %.1f sec.\n', brain_decimation.efield_elapsed_sec);

        % etc_tms_efield_surf returns vertices relative to coords. Replace them
        % with the original full-brain vertex numbers so sparse overlays still
        % address the subject's regular FreeSurfer brain surface.
        if(isfield(efield_sparse,'vertices'))
            brain_decimation.sparse_query_vertices=efield_sparse.vertices;
        end;
        efield_sparse.vertices=int32(brain_decimation.sparse_full_vertex_index(:)-1);

        efield_sparse.tms_coil_xfm=tms_coil_xfm_moved;
        efield_sparse.tms_coil_xfm_mm=tms_coil_xfm_mm;
        efield_sparse.coil_center=coil_center;
        efield_sparse.coil_center_native=coil_center;
        efield_sparse.coil_center_mni=coil_center_mni;
        efield_sparse.coil_orientation=coil_orientation;
        efield_sparse.coil_orientation_native=coil_orientation;
        efield_sparse.coil_orientation_mni=coil_orientation_mni;
        efield_sparse.coil_up=coil_up;
        efield_sparse.coil_up_native=coil_up;
        efield_sparse.coil_up_mni=coil_up_mni;
        efield_sparse.brainsight_AP_deg=brainsight_pose.AP_deg;
        efield_sparse.brainsight_Lat_deg=brainsight_pose.Lat_deg;
        efield_sparse.brainsight_Twist_deg=brainsight_pose.Twist_deg;
        efield_sparse.brainsight_AP_Lat_Twist_deg=[brainsight_pose.AP_deg brainsight_pose.Lat_deg brainsight_pose.Twist_deg];
        efield_sparse.coil_rotation_deg=target_coil_rotation_deg_this;
        efield_sparse.tms_coil_name=tms_coil_name;
        efield_sparse.basis_solve_elapsed_sec=basis_solve_elapsed_sec;
        efield_sparse.target_coord=target_coord_this;
        efield_sparse.target_coord_mni=target_coords_mni(target_idx,:);
        efield_sparse.target_source=target_source_this;
        efield_sparse.target_candidate_info=target_candidate_info;
        efield_sparse.target_index=target_idx;
        efield_sparse.target_count=n_targets;
        efield_sparse.output_stem=output_stem_target;
        efield_sparse.bem_obj=bem_obj;
        efield_sparse.scalp_bem_index=head_surf_idx;
        efield_sparse.decimation=brain_decimation;
        efield_sparse.flag_sparse=1;

        % for tms_nav window....
        etc_render_fsbrain.efield=efield_sparse;
        strcoil=strcoil_target;

        if(flag_nav)
            brain_decimation_render=brain_decimation;
            efield_render=etc_local_expand_efield_to_full(efield_sparse, brain_decimation_render, brain_decimation_render.full_nvertices, full_interpolation_smooth_step);
            render_view_state=etc_local_render_efield(efield_render, output_stem_target, overlay_threshold, render_view_state);

            % for tms_nav window....
            %enable efield entries 
            etc_render_fsbrain.app_tms_nav.EfieldShowCheckBox.Value=1;
            etc_render_fsbrain.app_tms_nav.EfieldDropDown.Value='Total';
            etc_render_fsbrain.app_tms_nav.efield.Etotal=efield_render.Etotal;
            etc_render_fsbrain.app_tms_nav.efield.Enormal=efield_render.Enormal;
            etc_render_fsbrain.app_tms_nav.efield.Etangent=efield_render.Etangent;
            etc_render_fsbrain.app_tms_nav.efield.vertices=efield_render.vertices;
            etc_local_set_nav_efield_status('g');
        end;

        if(flag_interpolate_to_full_brain)
            fprintf('Expanding sparse E-field to the full brain mesh with nearest fill + inverse_smooth step=%d...\n',full_interpolation_smooth_step);
            efield=etc_local_expand_efield_to_full(efield_sparse, brain_decimation, brain_decimation.full_nvertices, full_interpolation_smooth_step);
            efield_sparse.decimation=brain_decimation;
        else
            efield=efield_sparse;
        end;

        efield.decimation=brain_decimation;
        efield.coil_center_native=coil_center;
        efield.coil_center_mni=coil_center_mni;
        efield.coil_orientation_native=coil_orientation;
        efield.coil_orientation_mni=coil_orientation_mni;
        efield.coil_up_native=coil_up;
        efield.coil_up_mni=coil_up_mni;
        efield.brainsight_AP_deg=brainsight_pose.AP_deg;
        efield.brainsight_Lat_deg=brainsight_pose.Lat_deg;
        efield.brainsight_Twist_deg=brainsight_pose.Twist_deg;
        efield.brainsight_AP_Lat_Twist_deg=[brainsight_pose.AP_deg brainsight_pose.Lat_deg brainsight_pose.Twist_deg];
        efield.basis_solve_elapsed_sec=basis_solve_elapsed_sec;
        if(flag_interpolate_to_full_brain)
            efield_sparse.decimation=brain_decimation;
        end;

        if(flag_on_target)
            if(flag_interpolate_to_full_brain)
                brain_decimation_ontarget=brain_decimation;
                efield_ontarget=efield;
            else
                brain_decimation_ontarget=brain_decimation;
                efield_ontarget=etc_local_expand_efield_to_full(efield_sparse, brain_decimation_ontarget, brain_decimation_ontarget.full_nvertices, full_interpolation_smooth_step);
            end;

            on_target=etc_local_compute_on_target( ...
                efield_ontarget, brain_decimation_ontarget, target_coord_this, ...
                on_target_label_files, on_target_circle_diameter_mm, on_target_circle_hemi, ...
                subjects_dir, subject{subj_idx}, on_target_exclude_sulcal_vertices, ...
                on_target_sulc_threshold, on_target_sulc_keep_direction, ...
                on_target_empty_after_sulcal_mask_action, on_target_field_name, ...
                on_target_hotspot_percentile_thresholds, on_target_hotspot_hemi, ...
                flag_on_target_didt, on_target_didt_reference_A_per_us, ...
                on_target_didt_A_per_us, on_target_didt_activation_threshold_V_per_m, ...
                on_target_didt_hotspot_hemi);
            if(strcmpi(on_target.roi.method,'circle'))
                on_target=etc_local_write_on_target_circle_labels(on_target, brain_decimation_ontarget, root_dir, output_stem_target);
            end;
            if(flag_on_target_write_hotspot_label)
                on_target=etc_local_write_hotspot_efield_labels(on_target, efield_ontarget, brain_decimation_ontarget, root_dir, output_stem_target, on_target_hotspot_label_percentile);
            else
                on_target=etc_local_attach_hotspot_label_summary(on_target, on_target_hotspot_label_percentile);
            end;
            if(flag_on_target_didt && isfield(on_target,'didt'))
                on_target=etc_local_write_didt_dose_csv(on_target, root_dir, output_stem_target);
            end;
            efield.on_target=on_target;
            efield_sparse.on_target=on_target;
            brain_decimation.on_target=on_target.roi;

            fprintf('On-target ROI [%s]: n=%d, area=%1.1f mm^2, Etotal mean=%1.3f, median=%1.3f, max=%1.3f V/m\n', ...
                on_target.roi.method, on_target.roi.nvertices, on_target.roi.area_mm2, on_target.Etotal.mean, on_target.Etotal.median, on_target.Etotal.max);
            if(isfield(on_target.roi,'sulcal_mask'))
                if(isfield(on_target.roi.sulcal_mask,'flag_exclude'))
                    if(on_target.roi.sulcal_mask.flag_exclude)
                        if(isfield(on_target.roi.sulcal_mask,'flag_recentered') && on_target.roi.sulcal_mask.flag_recentered)
                            fprintf('Sulcal ROI mask: initial disk kept %d/%d vertices (%1.1f%%) using %s %1.3f; recentered to nearest gyral vertex %d (%1.2f mm from target), final n=%d (%1.1f%% of initial)\n', ...
                                on_target.roi.sulcal_mask.nvertices_after_pre_recenter, ...
                                on_target.roi.sulcal_mask.nvertices_before, on_target.roi.sulcal_mask.percent_after_pre_recenter, ...
                                on_target.roi.sulcal_mask.rule, on_target.roi.sulcal_mask.threshold, ...
                                on_target.roi.sulcal_mask.recenter_vertex_index1, on_target.roi.sulcal_mask.recenter_distance_to_target_mm, ...
                                on_target.roi.nvertices, on_target.roi.sulcal_mask.percent_after);
                        else
                            fprintf('Sulcal ROI mask: kept %d/%d vertices (%1.1f%%), removed %d using %s %1.3f\n', ...
                                on_target.roi.sulcal_mask.nvertices_after, on_target.roi.sulcal_mask.nvertices_before, ...
                                on_target.roi.sulcal_mask.percent_after, on_target.roi.sulcal_mask.nvertices_removed, ...
                                on_target.roi.sulcal_mask.rule, on_target.roi.sulcal_mask.threshold);
                        end;
                    end;
                end;
            end;
            fprintf('TANS on-target using %s hotspot percentiles [%1.1f..%1.1f]: mean=%1.2f%%\n', ...
                on_target.hotspot.field_name, max(on_target.hotspot.percentile_thresholds), min(on_target.hotspot.percentile_thresholds), 100.*on_target.TANS_OnTarget);
            if(isfield(on_target.hotspot,'label_percentile'))
                fprintf('Hotspot label p%1.1f: n=%d, area=%1.1f mm^2, target-overlap=%1.2f%% of hotspot\n', ...
                    on_target.hotspot.label_percentile, on_target.hotspot.label.nvertices, on_target.hotspot.label.area_mm2, 100.*on_target.hotspot.label.on_target);
            end;
            etc_local_print_didt_dose_on_target(on_target);
        end;

        target_summary_row=etc_local_target_summary_row(subject{subj_idx}, target_idx, output_stem_target, ...
            target_coord_this, target_coords_mni(target_idx,:), target_coil_rotation_deg_this, ...
            target_source_this, coil_center, coil_center_mni, coil_orientation, coil_orientation_mni, ...
            coil_up, coil_up_mni, brainsight_pose, efield, on_target, label_info);
        target_summary(end+1)=target_summary_row; %#ok<SAGROW>

        if(flag_stc)
            if(flag_interpolate_to_full_brain)
                efield_stc=efield;
                brain_decimation_stc=brain_decimation;
            elseif(~isempty(efield_ontarget))
                efield_stc=efield_ontarget;
                brain_decimation_stc=brain_decimation_ontarget;
            else
                fprintf('Expanding sparse E-field to the full brain mesh for native STC output with nearest fill + inverse_smooth step=%d...\n',full_interpolation_smooth_step);
                brain_decimation_stc=brain_decimation;
                efield_stc=etc_local_expand_efield_to_full(efield_sparse, brain_decimation_stc, brain_decimation_stc.full_nvertices, full_interpolation_smooth_step);
            end;
            stc_info=etc_local_write_efield_stc(efield_stc, brain_decimation_stc, root_dir, output_stem_target);
            efield.stc=stc_info;
            efield_sparse.stc=stc_info;
        end;

        native_nii_info=[];

        target_info=struct;
        target_info.subject=subject{subj_idx};
        target_info.subject_index=subj_idx;
        target_info.index=target_idx;
        target_info.count=n_targets;
        target_info.coord=target_coord_this;
        target_info.coord_mni=target_coords_mni(target_idx,:);
        target_info.all_coords=target_coords;
        target_info.all_coords_mni=target_coords_mni;
        target_info.output_stem=output_stem_target;
        target_info.target_source=target_source_this;
        target_info.target_candidate_info=target_candidate_info;
        target_info.coil_rotation_deg=target_coil_rotation_deg_this;
        target_info.flag_stc=flag_stc;
        target_info.stc=stc_info;
        target_info.tms_coil_xfm_goto=tms_coil_xfm_goto;
        target_info.tms_coil_xfm=tms_coil_xfm_moved;
        target_info.tms_coil_xfm_mm=tms_coil_xfm_mm;
        target_info.coil_center=coil_center;
        target_info.coil_center_native=coil_center;
        target_info.coil_center_mni=coil_center_mni;
        target_info.coil_orientation=coil_orientation;
        target_info.coil_orientation_native=coil_orientation;
        target_info.coil_orientation_mni=coil_orientation_mni;
        target_info.coil_up=coil_up;
        target_info.coil_up_native=coil_up;
        target_info.coil_up_mni=coil_up_mni;
        target_info.brainsight_AP_deg=brainsight_pose.AP_deg;
        target_info.brainsight_Lat_deg=brainsight_pose.Lat_deg;
        target_info.brainsight_Twist_deg=brainsight_pose.Twist_deg;
        target_info.brainsight_AP_Lat_Twist_deg=[brainsight_pose.AP_deg brainsight_pose.Lat_deg brainsight_pose.Twist_deg];
        target_info.native_nii=native_nii_info;

        result_cache_entry=struct;
        result_cache_entry.summary_row=target_summary_row;
        result_cache_entry.efield=efield;
        result_cache_entry.efield_sparse=efield_sparse;
        result_cache_entry.brain_decimation=brain_decimation;
        result_cache_entry.efield_ontarget=efield_ontarget;
        result_cache_entry.brain_decimation_ontarget=brain_decimation_ontarget;
        result_cache_entry.on_target=on_target;
        result_cache_entry.target_info=target_info;
        result_cache_entry.stc_info=stc_info;
        result_cache_entry.native_nii=native_nii_info;
        top_result_cache=etc_local_update_top_result_cache(top_result_cache, result_cache_entry, top_result_count);

        % Candidate rows are ordered by unique ROI position, then rotation.
        % Save and validate only after the final rotation of a position.
        is_last_rotation=(target_idx==n_targets) || ...
            (double(target_candidate_info.table.candidate_id(target_idx+1))~=candidate_id_this);

        if(flag_validate_two_basis_synthesis && candidate_id_this>0 && ...
                candidate_id_this<=validation_candidate_count && is_last_rotation)
            basis_validation_rows=[basis_validation_rows; etc_local_validate_basis_candidate( ...
                candidate_id_this, target_coord_this, candidate_coil_rotation_deg, ...
                tms_coil_xfm_goto, tms_coil_xfm, strcoil_base, ...
                bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, ...
                cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, ...
                RnumberE, ineighborE, EC, coords, tissue_to_plot, basis_cache{candidate_id_this})]; %#ok<AGROW>
        end;

        if(flag_write_search_checkpoint && candidate_id_this>0 && is_last_rotation && ...
                (mod(candidate_id_this,search_checkpoint_interval)==0 || target_idx==n_targets))
            next_target_idx=target_idx+1;
            save(search_checkpoint_file,'target_summary','top_result_cache','basis_cache', ...
                'next_target_idx','checkpoint_signature','-v7.3');
            fprintf('Saved two-basis checkpoint after ROI position %d/%d.\n', ...
                candidate_id_this,target_candidate_info.n_unique_coords);
        end;
    end;

    % Keep a separate, complete table for derived volumes.  The ranked/top
    % table below is only for navigation and reporting; it must never limit
    % the candidate/rotation values used for TANS volume export.
    all_target_summary_table=struct2table(target_summary);
    target_summary_table=sortrows(all_target_summary_table, ...
        {'TANS_OnTarget_percent','TANS_target_coverage_percent','Etotal_roi_mean_Vm','Etotal_roi_max_Vm'}, ...
        {'descend','descend','descend','descend'});
    n_top=min(top_result_count,etc_local_table_nrows(target_summary_table));
    top_target_table=target_summary_table(1:n_top,:);
    all_metrics_csv=fullfile(root_dir,sprintf('%s_all_target_metrics.csv',output_stem));
    writetable(target_summary_table,all_metrics_csv);
    if(flag_validate_two_basis_synthesis && ~isempty(basis_validation_rows))
        basis_validation_csv=fullfile(root_dir,sprintf('%s_two_basis_validation.csv',output_stem));
        writetable(struct2table(basis_validation_rows),basis_validation_csv);
        fprintf('Wrote two-basis validation results: %s\n',basis_validation_csv);
    end;
    tans_nii_info=[];
    if(flag_native_tans_nii)
        tans_nii_info=etc_local_write_tans_nii_outputs( ...
            subjects_dir, subject{subj_idx}, root_dir, output_stem, ...
            all_target_summary_table, bem_obj(head_surf_idx).vertex, nearest_chunk_size);
    end;
    tans_surface_info=etc_local_write_tans_surface_outputs( ...
        all_target_summary_table, scalp_vertex_full, root_dir, output_stem, nearest_chunk_size);
    top_metrics_csv=fullfile(root_dir,sprintf('%s_ranked_targets.csv',output_stem));
    writetable(top_target_table,top_metrics_csv);
    best_target=top_target_table(1,:);
    atlas_comparison=target_summary_table(etc_local_string_equal(target_summary_table.target_source,'atlas_mni305'),:);
    best_optimized=target_summary_table(etc_local_string_equal(target_summary_table.target_source,'dfconn_optimized_candidate'),:);
    if(etc_local_table_nrows(atlas_comparison)>0 && etc_local_table_nrows(best_optimized)>0)
        atlas_vs_best=etc_local_make_atlas_vs_best_table(atlas_comparison(1,:), best_optimized(1,:));
        atlas_vs_best_csv=fullfile(root_dir,sprintf('%s_atlas_vs_best_optimized.csv',output_stem));
        writetable(atlas_vs_best,atlas_vs_best_csv);
    else
        atlas_vs_best=table();
        atlas_vs_best_csv='';
    end;
    if(flag_nav)
        render_view_state=etc_local_show_best_nav_result( ...
            top_result_cache, best_target, render_view_state, overlay_threshold, ...
            full_interpolation_smooth_step);
    end;
    etc_local_save_top_result_cache(top_result_cache, root_dir, label_info, subjects_dir, flag_native_nii, flag_stc, target_coords, target_coords_mni, target_coord_mni, target_coil_rotation_deg, target_candidate_info, on_target_hotspot_label_percentile);
    fprintf('Optimal target native coordinate [%1.2f %1.2f %1.2f] mm; MNI coordinate [%1.2f %1.2f %1.2f] mm; coil rotation %1.2f deg; TANS on-target=%1.2f%%.\n', ...
        best_target.target_native_x, best_target.target_native_y, best_target.target_native_z, ...
        best_target.target_mni_x, best_target.target_mni_y, best_target.target_mni_z, ...
        best_target.coil_rotation_deg, best_target.TANS_OnTarget_percent);
    fprintf('Optimal coil scalp position native [%1.2f %1.2f %1.2f] mm; MNI [%1.2f %1.2f %1.2f] mm; native axis [%1.3f %1.3f %1.3f]; MNI axis [%1.3f %1.3f %1.3f].\n', ...
        best_target.coil_center_x, best_target.coil_center_y, best_target.coil_center_z, ...
        best_target.coil_center_mni_x, best_target.coil_center_mni_y, best_target.coil_center_mni_z, ...
        best_target.coil_axis_x, best_target.coil_axis_y, best_target.coil_axis_z, ...
        best_target.coil_axis_mni_x, best_target.coil_axis_mni_y, best_target.coil_axis_mni_z);
    fprintf('Brainsight pose: AP=%1.2f deg; Lat=%1.2f deg; Twist=%1.2f deg.\n', ...
        best_target.brainsight_AP_deg, best_target.brainsight_Lat_deg, best_target.brainsight_Twist_deg);
    fprintf('Wrote all tried target metrics: [%s]\n',all_metrics_csv);
    fprintf('Wrote top %d ranked target metrics: [%s]\n',n_top,top_metrics_csv);
    fprintf('Saved top %d e-field MAT files and hotspot labels.\n',n_top);
    if(~isempty(atlas_vs_best_csv))
        fprintf('Wrote atlas-vs-best-optimized comparison: [%s]\n',atlas_vs_best_csv);
    end;
    if(flag_write_search_checkpoint && exist(search_checkpoint_file,'file')==2)
        delete(search_checkpoint_file);
        fprintf('Removed completed search checkpoint [%s].\n',search_checkpoint_file);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view_state=etc_local_show_best_nav_result(top_result_cache, best_target, view_state, overlay_threshold, smooth_step)

global etc_render_fsbrain

best_target_idx=double(best_target.target_index(1));
best_cache_idx=[];
for cache_idx=1:numel(top_result_cache)
    if(double(top_result_cache{cache_idx}.target_info.index)==best_target_idx)
        best_cache_idx=cache_idx;
        break;
    end;
end;
if(isempty(best_cache_idx))
    warning('Cannot find ranked target %d in the navigation result cache.',best_target_idx);
    return;
end;

best_entry=top_result_cache{best_cache_idx};
best_target_info=best_entry.target_info;
etc_local_apply_coil_to_nav(best_target_info.tms_coil_xfm);
etc_local_update_nav_target_coordinate(best_target_info.coord);

efield_best=best_entry.efield;
if(~isfield(efield_best,'flag_sparse') || efield_best.flag_sparse~=0)
    efield_best=etc_local_expand_efield_to_full( ...
        best_entry.efield_sparse, best_entry.brain_decimation, ...
        best_entry.brain_decimation.full_nvertices, smooth_step);
end;
efield_best.decimation=best_entry.brain_decimation;
etc_render_fsbrain.efield=efield_best;

try
    if(isfield(etc_render_fsbrain,'app_tms_nav') && ...
            ~isempty(etc_render_fsbrain.app_tms_nav) && ...
            isvalid(etc_render_fsbrain.app_tms_nav))
        etc_render_fsbrain.app_tms_nav.EfieldShowCheckBox.Value=1;
        etc_render_fsbrain.app_tms_nav.EfieldDropDown.Value='Total';
        etc_render_fsbrain.app_tms_nav.efield.Etotal=efield_best.Etotal;
        etc_render_fsbrain.app_tms_nav.efield.Enormal=efield_best.Enormal;
        etc_render_fsbrain.app_tms_nav.efield.Etangent=efield_best.Etangent;
        etc_render_fsbrain.app_tms_nav.efield.vertices=efield_best.vertices;
        etc_local_set_nav_efield_status('g');
    end;
catch
end;

view_state=etc_local_render_efield( ...
    efield_best, char(best_target_info.output_stem), overlay_threshold, view_state);
fprintf('Navigation GUI updated to optimal target %d [%s].\n', ...
    best_target_idx, char(best_target_info.output_stem));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_update_nav_target_coordinate(target_coord)

global etc_render_fsbrain

try
    if(~isfield(etc_render_fsbrain,'app_tms_nav') || ...
            isempty(etc_render_fsbrain.app_tms_nav) || ...
            ~isvalid(etc_render_fsbrain.app_tms_nav))
        return;
    end;

    target_coord=double(target_coord(:).');
    if(numel(target_coord)~=3)
        return;
    end;

    % Match etc_tms_init_nav/ etc_render_fsbrain_handle: identify the
    % closest surface vertex and derive CRS, MNI305, and Scanner coordinates
    % from the current etc_render_fsbrain volume geometry.
    vertex_coord=double(etc_render_fsbrain.vertex_coords);
    dist2=sum((vertex_coord-repmat(target_coord,size(vertex_coord,1),1)).^2,2);
    [~,target_vertex_index]=min(dist2);

    click_vertex_vox=inv(etc_render_fsbrain.vol.tkrvox2ras)*[target_coord 1].';
    click_vertex_vox=click_vertex_vox(1:3).';
    click_vertex_point_mni=etc_render_fsbrain.talxfm*etc_render_fsbrain.vol_pre_xfm* ...
        etc_render_fsbrain.vol.vox2ras*[click_vertex_vox 1].';
    click_vertex_point_mni=click_vertex_point_mni(1:3).';
    scanner_coord=etc_render_fsbrain.vol.vox2ras*[click_vertex_vox 1].';
    scanner_coord=scanner_coord(1:3).';

    % Keep the renderer's click state synchronized as well as the app fields.
    etc_render_fsbrain.click_coord=target_coord;
    etc_render_fsbrain.click_vertex_vox=click_vertex_vox;
    etc_render_fsbrain.click_vertex_vox_round=round(click_vertex_vox);
    etc_render_fsbrain.click_vertex_point_tal=click_vertex_point_mni;
    etc_render_fsbrain.click_vertex_point_round_tal=click_vertex_point_mni;

    app=etc_render_fsbrain.app_tms_nav;
    etc_render_fsbrain_tms_nav_notify(app,struct('Source',app.vertexindexEditField),target_vertex_index);
    etc_render_fsbrain_tms_nav_notify(app,struct('Source',app.XYZEditField),target_coord);
    etc_render_fsbrain_tms_nav_notify(app,struct('Source',app.CRSEditField),click_vertex_vox);
    etc_render_fsbrain_tms_nav_notify(app,struct('Source',app.MNIEditField),click_vertex_point_mni);
    etc_render_fsbrain_tms_nav_notify(app,struct('Source',app.ScannerEditField),scanner_coord);
catch ME
    fprintf('Navigation target coordinate update failed: %s\n',ME.message);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_set_nav_efield_status(color)

global etc_render_fsbrain

try
    if(isfield(etc_render_fsbrain,'app_tms_nav') && ...
            ~isempty(etc_render_fsbrain.app_tms_nav) && ...
            isvalid(etc_render_fsbrain.app_tms_nav))
        app=etc_render_fsbrain.app_tms_nav;
        etc_render_fsbrain_tms_nav_notify( ...
            app,struct('Source',app.EfieldCalclLamp),color);
        drawnow;
    end;
catch
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_info=etc_local_create_dfconn_sup_mPFC_labels(subject, tms_dir, volume_file, decimation, config, batch_log_dir)

if(~exist(tms_dir,'dir'))
    mkdir(tms_dir);
end;

subject_label_files=etc_local_ensure_subject_sup_mPFC_labels(subject, tms_dir, config, batch_log_dir);

label_info=struct;
label_info.subject=subject;
label_info.volume_file=volume_file;
label_info.volume_threshold=config.volume_threshold;
label_info.fsaverage_label_files={config.fsaverage_lh_label,config.fsaverage_rh_label};
label_info.subject_sup_mPFC_label_files=subject_label_files;
label_info.output_prefix=fullfile(tms_dir,config.label_prefix);
label_info.hemi=struct([]);
label_info.global_vertex_index1=[];
label_info.global_vertex_index0=int32([]);
label_info.dfconn_value=[];
label_info.label_files={};

for hemi_idx=1:numel(decimation.hemi)
    hemi_name=lower(decimation.hemi(hemi_idx).name);
    if(~strcmp(hemi_name,'lh') && ~strcmp(hemi_name,'rh'))
        continue;
    end;

    surface_value=etc_local_sample_volume_to_surface(subject, volume_file, hemi_name, config);
    full_count=decimation.hemi(hemi_idx).full_vertex_count;
    if(numel(surface_value)~=full_count)
        error('mri_vol2surf returned %d %s vertices, expected %d.',numel(surface_value),hemi_name,full_count);
    end;

    base_label_file=subject_label_files.(hemi_name);
    base_local_idx1=etc_local_read_label_local_indices(base_label_file, full_count);
    base_mask=false(full_count,1);
    base_mask(base_local_idx1)=true;

    keep_local_idx1=find(base_mask & isfinite(surface_value(:)) & surface_value(:)>config.volume_threshold);
    full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
    keep_global_idx1=keep_local_idx1(:)+full_offset;
    keep_value=surface_value(keep_local_idx1);

    label_file=sprintf('%s-%s.label',label_info.output_prefix,hemi_name);
    vertex_coord=decimation.full_vertex_coords(keep_global_idx1,:);
    inverse_write_label(keep_local_idx1(:)-1, vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3), keep_value(:), label_file);

    label_info.hemi(hemi_idx).name=hemi_name;
    label_info.hemi(hemi_idx).label_file=label_file;
    label_info.hemi(hemi_idx).sampled_value_nvertices=numel(surface_value);
    label_info.hemi(hemi_idx).base_label_file=base_label_file;
    label_info.hemi(hemi_idx).base_label_nvertices=numel(base_local_idx1);
    label_info.hemi(hemi_idx).dfconn_gt_threshold_nvertices=sum(isfinite(surface_value(:)) & surface_value(:)>config.volume_threshold);
    label_info.hemi(hemi_idx).nvertices=numel(keep_local_idx1);
    label_info.hemi(hemi_idx).percent_of_base_label=100.*numel(keep_local_idx1)./max(1,numel(base_local_idx1));
    label_info.hemi(hemi_idx).local_vertex_index1=keep_local_idx1(:);
    label_info.hemi(hemi_idx).local_vertex_index0=int32(keep_local_idx1(:)-1);
    label_info.hemi(hemi_idx).global_vertex_index1=keep_global_idx1(:);
    label_info.hemi(hemi_idx).global_vertex_index0=int32(keep_global_idx1(:)-1);
    label_info.hemi(hemi_idx).dfconn_value=keep_value(:);
    label_info.hemi(hemi_idx).dfconn_summary=etc_local_value_summary(keep_value(:));

    label_info.global_vertex_index1=cat(1,label_info.global_vertex_index1(:),keep_global_idx1(:));
    label_info.global_vertex_index0=int32(label_info.global_vertex_index1(:)-1);
    label_info.dfconn_value=cat(1,label_info.dfconn_value(:),keep_value(:));
    label_info.label_files{end+1,1}=label_file;
end;

label_info.total_nvertices=numel(label_info.global_vertex_index1);
label_info.dfconn_summary=etc_local_value_summary(label_info.dfconn_value);
label_info.combined_label_file=sprintf('%s_both-lh.label',label_info.output_prefix);
combined_coord=decimation.full_vertex_coords(label_info.global_vertex_index1,:);
inverse_write_label(label_info.global_vertex_index0(:), combined_coord(:,1), combined_coord(:,2), combined_coord(:,3), label_info.dfconn_value(:), label_info.combined_label_file);
label_info.label_files{end+1,1}=label_info.combined_label_file;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function native_nii_info=etc_local_write_native_nii_outputs(subjects_dir, subject, efield_full, label_files, output_dir, output_stem)

if(exist('MRIwrite','file')~=2)
    error('Cannot find MRIwrite.m on the MATLAB path; needed for native NIfTI export.');
end;

if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end;

template_file=fullfile(subjects_dir,subject,'mri','orig.mgz');
if(exist(template_file,'file')~=2)
    error('Cannot find native volume template [%s].',template_file);
end;

native_nii_info=struct;
native_nii_info.template_file=template_file;
native_nii_info.output_dir=output_dir;
native_nii_info.output_stem=output_stem;
native_nii_info.efield=etc_local_write_native_efield_nii(template_file, subject, efield_full, output_dir, output_stem);
native_nii_info.roi=etc_local_write_native_roi_nii(template_file, subject, label_files, output_dir, output_stem);

fprintf('Wrote native NIfTI outputs for %s: ROI=%s, E-field=%s.\n', ...
    subject, native_nii_info.roi.both_file, native_nii_info.efield.Etotal_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function efield_nii_info=etc_local_write_native_efield_nii(template_file, subject, efield_full, output_dir, output_stem)

scalar_names={'Etotal','Enormal','Etangent'};
efield_nii_info=struct;
efield_nii_info.output_dir=output_dir;
efield_nii_info.output_stem=output_stem;

for scalar_idx=1:length(scalar_names)
    scalar_name=scalar_names{scalar_idx};
    value=etc_local_get_efield_stc_value(efield_full, scalar_name);
    decimation=efield_full.decimation;
    if(numel(value)==decimation.full_nvertices)
        full_value=double(value(:));
    elseif(numel(value)==decimation.sparse_nvertices)
        full_value=nan(decimation.full_nvertices,1);
        full_value(decimation.sparse_full_vertex_index(:))=double(value(:));
    else
        error('%s length %d does not match full_nvertices %d or sparse_nvertices %d.', ...
            scalar_name, numel(value), decimation.full_nvertices, decimation.sparse_nvertices);
    end;
    full_value(~isfinite(full_value))=0;

    out_file=fullfile(output_dir,sprintf('%s_native_%s.nii.gz',output_stem,scalar_name));
    [lh_value, rh_value]=etc_local_split_hemi_values(decimation, full_value);
    etc_local_write_native_surface_volume(subjects_dir_from_template(template_file), subject, template_file, lh_value, rh_value, out_file);
    efield_nii_info.(scalar_name)=out_file;
    fprintf('Wrote native NIfTI E-field [%s]: %s\n',scalar_name,out_file);
end;

efield_nii_info.Etotal_file=efield_nii_info.Etotal;
efield_nii_info.Enormal_file=efield_nii_info.Enormal;
efield_nii_info.Etangent_file=efield_nii_info.Etangent;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi_nii_info=etc_local_write_native_roi_nii(template_file, subject, label_files, output_dir, output_stem)

if(isstruct(label_files))
    if(isfield(label_files,'lh') && isfield(label_files,'rh'))
        label_files={label_files.lh,label_files.rh};
    else
        error('Expected label_files struct with lh and rh fields.');
    end;
end;
if(ischar(label_files))
    label_files={label_files};
end;
if(numel(label_files)<2)
    error('Expected lh and rh label files for native ROI export.');
end;

roi_nii_info=struct;
roi_nii_info.output_dir=output_dir;
roi_nii_info.output_stem=output_stem;
roi_nii_info.template_file=template_file;

subjects_dir=fileparts(fileparts(fileparts(template_file)));
hemi={'lh','rh'};
roi_lh_value=[];
roi_rh_value=[];

for hemi_idx=1:numel(hemi)
    hemi_name=hemi{hemi_idx};
    label_file=label_files{hemi_idx};
    label_data=etc_local_read_label_file(subject, label_file);
    if(size(label_data,2)<4)
        error('Label file [%s] does not contain vertex coordinates.',label_file);
    end;

    hemi_full_count=etc_local_surface_vertex_count(subjects_dir, subject, hemi_name);
    hemi_full_values=zeros(hemi_full_count,1,'double');
    hemi_local_idx=label_data(:,1)+1;
    hemi_local_idx=hemi_local_idx(isfinite(hemi_local_idx));
    hemi_local_idx=hemi_local_idx(hemi_local_idx>=1 & hemi_local_idx<=hemi_full_count);
    hemi_full_values(hemi_local_idx)=1;
    hemi_out_file=fullfile(output_dir,sprintf('%s_native_roi-%s.nii.gz',output_stem,hemi_name));
    if(strcmp(hemi_name,'lh'))
        roi_lh_value=hemi_full_values;
    else
        roi_rh_value=hemi_full_values;
    end;
    roi_nii_info.(hemi_name)=hemi_out_file;
    fprintf('Wrote native ROI NIfTI [%s]: %s\n',hemi_name,hemi_out_file);
end;
roi_both_file=fullfile(output_dir,sprintf('%s_native_roi-both.nii.gz',output_stem));
etc_local_write_native_surface_volume(subjects_dir, subject, template_file, roi_lh_value, roi_rh_value, roi_both_file);
roi_nii_info.both=roi_both_file;
roi_nii_info.both_file=roi_both_file;
fprintf('Wrote native ROI NIfTI [both]: %s\n',roi_both_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vol_out=etc_local_surface_points_to_native_volume(template_file, point_coord, point_value)

if(nargin<3 || isempty(point_value))
    point_value=ones(size(point_coord,1),1);
end;
point_coord=double(point_coord);
point_value=double(point_value(:));
if(size(point_coord,1)~=numel(point_value))
    error('point_coord and point_value must have matching row counts.');
end;

vol_out=MRIread(template_file);
if(~isfield(vol_out,'tkrvox2ras') || isempty(vol_out.tkrvox2ras))
    error('Native volume template is missing voxel geometry.');
end;

vol_size=size(vol_out.vol);
if(numel(vol_size)<3)
    error('Native volume template must be at least 3D.');
end;
vol_size=vol_size(1:3);
vol=zeros(vol_size,'single');

if(~isempty(point_coord))
    % Surface vertex coordinates are in FreeSurfer/tkr RAS.  They must be
    % mapped through tkrvox2ras, not vox2ras (which is scanner RAS).  Using
    % vox2ras produces plausible-looking but anatomically displaced points.
    vox=inv(vol_out.tkrvox2ras)*[point_coord ones(size(point_coord,1),1)].';
    vox=round(vox(1:3,:)).'+1;
    valid=isfinite(point_value) & all(isfinite(vox),2) & vox(:,1)>=1 & vox(:,1)<=vol_size(1) & vox(:,2)>=1 & vox(:,2)<=vol_size(2) & vox(:,3)>=1 & vox(:,3)<=vol_size(3);
    if(any(valid))
        lin=sub2ind(vol_size, vox(valid,1), vox(valid,2), vox(valid,3));
        [u_lin,~,ic]=unique(lin);
        agg_value=accumarray(ic, point_value(valid), [], @max);
        vol(u_lin)=single(agg_value);
    end;
end;

vol_out.vol=vol;
vol_out.nframes=1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lh_value, rh_value]=etc_local_split_hemi_values(decimation, full_value)

lh_value=zeros(decimation.hemi(1).full_vertex_count,1,'double');
rh_value=zeros(decimation.hemi(2).full_vertex_count,1,'double');

lh_off=decimation.hemi(1).full_vertex_offset;
lh_cnt=decimation.hemi(1).full_vertex_count;
rh_off=decimation.hemi(2).full_vertex_offset;
rh_cnt=decimation.hemi(2).full_vertex_count;

lh_value(:)=full_value((lh_off+1):(lh_off+lh_cnt));
rh_value(:)=full_value((rh_off+1):(rh_off+rh_cnt));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subjects_dir=subjects_dir_from_template(template_file)

subjects_dir=fileparts(fileparts(fileparts(template_file)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_write_native_surface_volume(subjects_dir, subject, template_file, lh_value, rh_value, out_file)

if(exist('write_curv','file')~=2)
    error('Cannot find write_curv on the MATLAB path.');
end;

lh_white_file=fullfile(subjects_dir,subject,'surf','lh.white');
rh_white_file=fullfile(subjects_dir,subject,'surf','rh.white');
lh_pial_file=fullfile(subjects_dir,subject,'surf','lh.pial');
rh_pial_file=fullfile(subjects_dir,subject,'surf','rh.pial');
if(exist(lh_white_file,'file')~=2 || exist(rh_white_file,'file')~=2 || ...
        exist(lh_pial_file,'file')~=2 || exist(rh_pial_file,'file')~=2)
    error('Cannot find white/pial subject surfaces for %s.',subject);
end;

lh_curv_file=[tempname '.curv'];
rh_curv_file=[tempname '.curv'];
cleanup_lh=onCleanup(@() etc_local_delete_if_exists(lh_curv_file)); %#ok<NASGU>
cleanup_rh=onCleanup(@() etc_local_delete_if_exists(rh_curv_file)); %#ok<NASGU>

write_curv(lh_curv_file, double(lh_value(:)), etc_local_surface_face_count(subjects_dir, subject, 'lh'));
write_curv(rh_curv_file, double(rh_value(:)), etc_local_surface_face_count(subjects_dir, subject, 'rh'));

surf2vol_bin=etc_local_find_executable('mri_surf2vol');
cmd=sprintf(['%s --sd %s --subject %s --template %s ' ...
    '--so %s %s --so %s %s --so %s %s --so %s %s --o %s'], ...
    surf2vol_bin, subjects_dir, subject, template_file, ...
    lh_white_file, lh_curv_file, lh_pial_file, lh_curv_file, ...
    rh_white_file, rh_curv_file, rh_pial_file, rh_curv_file, out_file);
[status,cmdout]=system(cmd);
if(status~=0)
    error('mri_surf2vol failed while writing [%s].\n%s',out_file,cmdout);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nfaces=etc_local_surface_face_count(subjects_dir, subject, hemi_name)

surf_file=fullfile(subjects_dir,subject,'surf',sprintf('%s.orig',hemi_name));
[~,faces]=read_surf(surf_file);
nfaces=size(faces,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nvertices=etc_local_surface_vertex_count(subjects_dir, subject, hemi_name)

% Label vertex indices are native FreeSurfer surface indices, so derive the
% expected overlay length directly from the corresponding subject surface.
surf_file=fullfile(subjects_dir,subject,'surf',sprintf('%s.orig',hemi_name));
if(exist(surf_file,'file')~=2)
    error('Cannot find subject surface file [%s].',surf_file);
end;
[vertices,~]=read_surf(surf_file);
nvertices=size(vertices,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function executable=etc_local_find_executable(name)

[status,cmdout]=system(sprintf('command -v %s',name));
cmdout=strtrim(cmdout);
if(status==0 && ~isempty(cmdout))
    executable=cmdout;
    return;
end;
error('Cannot find executable [%s].',name);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_write_native_volume_nifti(vol_out, template_file, out_file)

MRIwrite(vol_out,out_file);
fprintf('writing to %s (native header from %s)...\n', out_file, template_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tans_nii_info=etc_local_write_tans_nii_outputs(subjects_dir, subject, output_dir, output_stem, all_target_summary_table, scalp_vertex_coords, nearest_chunk_size)
if(exist('MRIread','file')~=2 || exist('MRIwrite','file')~=2)
    error('MRIread/MRIwrite are required for native TANS volume export.');
end;
template_file=fullfile(subjects_dir,subject,'mri','orig.mgz');
if(exist(template_file,'file')~=2)
    error('Cannot find native volume template [%s].',template_file);
end;
required_fields={'coil_center_x','coil_center_y','coil_center_z','coil_rotation_deg','TANS_OnTarget_percent'};
for field_idx=1:length(required_fields)
    if(~ismember(required_fields{field_idx},all_target_summary_table.Properties.VariableNames))
        error('TANS summary table is missing [%s].',required_fields{field_idx});
    end;
end;
target_coord=[double(all_target_summary_table.coil_center_x(:)) double(all_target_summary_table.coil_center_y(:)) double(all_target_summary_table.coil_center_z(:))];
rotation=double(all_target_summary_table.coil_rotation_deg(:));
tans=double(all_target_summary_table.TANS_OnTarget_percent(:));
valid=isfinite(target_coord(:,1)) & isfinite(target_coord(:,2)) & isfinite(target_coord(:,3)) & isfinite(rotation) & isfinite(tans);
if(~any(valid)) error('No finite target/TANS rows are available for native TANS volume export.'); end;
target_coord=target_coord(valid,:);
rotation=rotation(valid);
tans=tans(valid);
fprintf('Exporting native TANS volumes from %d/%d candidate/rotation rows.\n', ...
    numel(tans),height(all_target_summary_table));

if(isempty(scalp_vertex_coords) || size(scalp_vertex_coords,2)~=3)
    error('Decimated scalp BEM vertices are required for native TANS volume export.');
end;
scalp_vertex=etc_local_nearest_vertex_index(scalp_vertex_coords,target_coord,nearest_chunk_size);
[unique_vertex,~,target_group]=unique(scalp_vertex(:));
unique_coord=scalp_vertex_coords(unique_vertex,:);
rotation_values=unique(rotation(:)).';
score_by_rotation=nan(length(unique_vertex),length(rotation_values));
for row_idx=1:length(tans)
    rotation_idx=find(abs(rotation_values-rotation(row_idx))<=eps(max(1,abs(rotation(row_idx)))),1);
    current_value=score_by_rotation(target_group(row_idx),rotation_idx);
    if(isnan(current_value) || tans(row_idx)>current_value)
        score_by_rotation(target_group(row_idx),rotation_idx)=tans(row_idx);
    end;
end;
max_tans=nan(size(score_by_rotation,1),1);
mean_tans=nan(size(score_by_rotation,1),1);
best_rotation=nan(size(score_by_rotation,1),1);
for target_idx=1:size(score_by_rotation,1)
    value=score_by_rotation(target_idx,:);
    finite_idx=isfinite(value);
    if(any(finite_idx))
        max_tans(target_idx)=max(value(finite_idx));
        mean_tans(target_idx)=mean(value(finite_idx));
        candidate_rotation=rotation_values(finite_idx);
        candidate_value=value(finite_idx);
        [~,best_idx]=max(candidate_value);
        best_rotation(target_idx)=candidate_rotation(best_idx);
    end;
end;

template_vol=MRIread(template_file);
volume_size=size(template_vol.vol);
volume_size=volume_size(1:3);
all_rotation_vol=nan([volume_size length(rotation_values)],'single');
for rotation_idx=1:length(rotation_values)
    vol=etc_local_surface_points_to_native_volume(template_file,unique_coord,score_by_rotation(:,rotation_idx));
    all_rotation_vol(:,:,:,rotation_idx)=single(vol.vol);
end;

tans_nii_info=struct;
tans_nii_info.template_file=template_file;
tans_nii_info.output_dir=output_dir;
tans_nii_info.output_stem=output_stem;
tans_nii_info.rotation_deg=rotation_values;
tans_nii_info.n_target_locations=length(unique_vertex);
tans_nii_info.rotation_frame_csv=fullfile(output_dir,sprintf('%s_native_TANS_rotation_frames.csv',output_stem));
rotation_frame_table=table((1:length(rotation_values)).',rotation_values(:),'VariableNames',{'frame','coil_rotation_deg'});
writetable(rotation_frame_table,tans_nii_info.rotation_frame_csv);

% Preserve the exact complete input used to populate the volumes.  This is
% especially useful when several rotations share one target location.
candidate_value_table=table(target_coord(:,1),target_coord(:,2),target_coord(:,3), ...
    rotation,tans,scalp_vertex, ...
    'VariableNames',{'target_native_x','target_native_y','target_native_z', ...
    'coil_rotation_deg','TANS_OnTarget_percent','nearest_scalp_vertex_index'});
tans_nii_info.candidate_value_csv=fullfile(output_dir,sprintf('%s_native_TANS_candidate_values.csv',output_stem));
writetable(candidate_value_table,tans_nii_info.candidate_value_csv);

all_vol=template_vol;
all_vol.vol=all_rotation_vol;
all_vol.nframes=length(rotation_values);
all_vol.volsize=size(all_rotation_vol);
tans_nii_info.all_rotations=fullfile(output_dir,sprintf('%s_native_TANS_all_rotations.nii.gz',output_stem));
etc_local_write_native_volume_nifti(all_vol,template_file,tans_nii_info.all_rotations);
tans_nii_info.max=etc_local_write_tans_scalar_volume(template_file,unique_coord,max_tans,output_dir,output_stem,'max');
tans_nii_info.mean=etc_local_write_tans_scalar_volume(template_file,unique_coord,mean_tans,output_dir,output_stem,'mean');
tans_nii_info.best_rotation=etc_local_write_tans_scalar_volume(template_file,unique_coord,best_rotation,output_dir,output_stem,'best_rotation');
fprintf('Wrote native TANS volumes for %s: %s\n',subject,tans_nii_info.all_rotations);
end

function out_file=etc_local_write_tans_scalar_volume(template_file,target_coord,value,output_dir,output_stem,label)
vol=etc_local_surface_points_to_native_volume(template_file,target_coord,value);
out_file=fullfile(output_dir,sprintf('%s_native_TANS_%s.nii.gz',output_stem,label));
etc_local_write_native_volume_nifti(vol,template_file,out_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info=etc_local_write_tans_surface_outputs(summary_table, scalp_vertex_coords, output_dir, output_stem, nearest_chunk_size)

if(exist('write_wfile','file')~=2)
    warning('write_wfile is unavailable; skipping FreeSurfer scalp TANS overlays.');
    info=struct;
    return;
end;
required={'coil_center_x','coil_center_y','coil_center_z','coil_rotation_deg','TANS_OnTarget_percent'};
for idx=1:numel(required)
    if(~ismember(required{idx},summary_table.Properties.VariableNames))
        error('TANS scalp overlay table is missing [%s].',required{idx});
    end;
end;
coil_center=[double(summary_table.coil_center_x(:)) double(summary_table.coil_center_y(:)) double(summary_table.coil_center_z(:))];
rotation=double(summary_table.coil_rotation_deg(:));
tans=double(summary_table.TANS_OnTarget_percent(:));
valid=all(isfinite(coil_center),2) & isfinite(rotation) & isfinite(tans);
coil_center=coil_center(valid,:);
rotation=rotation(valid);
tans=tans(valid);
scalp_vertex=etc_local_nearest_vertex_index(double(scalp_vertex_coords),coil_center,nearest_chunk_size);
rotation_values=unique(rotation(:)).';
score=nan(size(scalp_vertex_coords,1),numel(rotation_values));
for idx=1:numel(tans)
    rotation_idx=find(abs(rotation_values-rotation(idx))<=eps(max(1,abs(rotation(idx)))),1);
    current=score(scalp_vertex(idx),rotation_idx);
    if(isnan(current) || tans(idx)>current)
        score(scalp_vertex(idx),rotation_idx)=tans(idx);
    end;
end;

max_value=nan(size(score,1),1);
mean_value=nan(size(score,1),1);
best_rotation=nan(size(score,1),1);
for vertex_idx=1:size(score,1)
    finite_idx=isfinite(score(vertex_idx,:));
    if(any(finite_idx))
        value=score(vertex_idx,finite_idx);
        rotations=rotation_values(finite_idx);
        max_value(vertex_idx)=max(value);
        mean_value(vertex_idx)=mean(value);
        [~,best_idx]=max(value);
        best_rotation(vertex_idx)=rotations(best_idx);
    end;
end;

info=struct;
info.max=fullfile(output_dir,sprintf('%s_native_TANS_outer_skin_max.w',output_stem));
info.mean=fullfile(output_dir,sprintf('%s_native_TANS_outer_skin_mean.w',output_stem));
info.best_rotation=fullfile(output_dir,sprintf('%s_native_TANS_outer_skin_best_rotation.w',output_stem));
etc_local_write_w_overlay(info.max,max_value);
etc_local_write_w_overlay(info.mean,mean_value);
etc_local_write_w_overlay(info.best_rotation,best_rotation);
info.rotation_files=cell(numel(rotation_values),1);
for rotation_idx=1:numel(rotation_values)
    info.rotation_files{rotation_idx}=fullfile(output_dir,sprintf('%s_native_TANS_outer_skin_rotation_%+05.1f.w',output_stem,rotation_values(rotation_idx)));
    etc_local_write_w_overlay(info.rotation_files{rotation_idx},score(:,rotation_idx));
end;
mapping_file=fullfile(output_dir,sprintf('%s_native_TANS_outer_skin_rotation_frames.csv',output_stem));
mapping_table=table((1:numel(rotation_values)).',rotation_values(:),'VariableNames',{'frame','coil_rotation_deg'});
writetable(mapping_table,mapping_file);
info.rotation_frame_csv=mapping_file;
fprintf('Wrote scalp TANS overlays on %d scalp vertices.\n',sum(any(isfinite(score),2)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_write_w_overlay(file_name,value)

value=double(value(:));
finite_idx=find(isfinite(value) & value~=0);
write_wfile(file_name,value(finite_idx),finite_idx-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_files=etc_local_ensure_subject_sup_mPFC_labels(subject, tms_dir, config, batch_log_dir)

setenv('SUBJECTS_DIR',config.subjects_dir);
label_files=struct;
label_files.lh=fullfile(tms_dir,sprintf('sup_mPFC_hippo_ent_fsaverage_to_%s-lh.label',subject));
label_files.rh=fullfile(tms_dir,sprintf('sup_mPFC_hippo_ent_fsaverage_to_%s-rh.label',subject));

if(exist(label_files.lh,'file')==2 && exist(label_files.rh,'file')==2)
    return;
end;

src_label=struct('lh',config.fsaverage_lh_label,'rh',config.fsaverage_rh_label);
hemis={'lh','rh'};
for idx=1:numel(hemis)
    hemi_name=hemis{idx};
    if(exist(label_files.(hemi_name),'file')==2)
        continue;
    end;
    log_file=fullfile(batch_log_dir,sprintf('%s_%s_sup_mPFC_hippo_ent_label2label_dfconn.log',subject,hemi_name));
    % Do not use POSIX shell redirection here. The server runs tcsh, where
    % "> logfile 2>&1" produces "Ambiguous output redirect". Capture the
    % command output through MATLAB instead, then write it to the batch log.
    cmd=sprintf(['mri_label2label --srcsubject fsaverage ' ...
        '--srclabel %s --trgsubject %s --trglabel %s --regmethod surface --hemi %s'], ...
        etc_local_shell_quote(src_label.(hemi_name)), ...
        etc_local_shell_quote(subject), ...
        etc_local_shell_quote(label_files.(hemi_name)), ...
        etc_local_shell_quote(hemi_name));
    [status,cmdout]=system(cmd);
    fid=fopen(log_file,'w');
    if(fid<0)
        warning('Cannot write mri_label2label log [%s].',log_file);
    else
        fprintf(fid,'%s\n',cmdout);
        fclose(fid);
    end;
    if(status~=0)
        error('mri_label2label failed for %s %s. See %s.',subject,hemi_name,log_file);
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surface_value=etc_local_sample_volume_to_surface(subject, volume_file, hemi_name, config)

setenv('SUBJECTS_DIR',config.subjects_dir);
out_file=sprintf('%s_%s_%s.mgh',tempname,subject,hemi_name);
cleanup=onCleanup(@() etc_local_delete_if_exists(out_file)); %#ok<NASGU>
cmd=sprintf(['mri_vol2surf --src %s --regheader %s --hemi %s ' ...
    '--surf %s --projfrac %g --interp %s --o %s'], ...
    etc_local_shell_quote(volume_file), ...
    etc_local_shell_quote(subject), ...
    etc_local_shell_quote(hemi_name), ...
    etc_local_shell_quote(config.surface_name), ...
    config.projfrac, ...
    etc_local_shell_quote(config.vol2surf_interp), ...
    etc_local_shell_quote(out_file));
[status,cmdout]=system(cmd);
if(status~=0)
    error('mri_vol2surf failed for %s %s: %s',subject,hemi_name,cmdout);
end;
m=MRIread(out_file);
surface_value=double(m.vol(:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target_coords,target_coords_mni,target_coil_rotation_deg,candidate_info,native_to_mni_xfm]=etc_local_build_tans_scalp_candidates(subjects_dir, subject, label_info, decimation, scalp_vertex, scalp_face, search_radius_mm, grid_spacing_mm, rotation_deg)

roi_idx=label_info.global_vertex_index1(:);
roi_coord=double(decimation.full_vertex_coords(roi_idx,:));
roi_value=double(label_info.dfconn_value(:));
weight=max(roi_value,0);
if(sum(weight)<=0 || any(~isfinite(weight)))
    weight=ones(size(weight));
end;
target_centroid=sum(roi_coord.*weight,1)./sum(weight);

scalp_vertex=double(scalp_vertex);
scalp_face=double(scalp_face);
scalp_normal=etc_local_surface_vertex_normals(scalp_vertex,scalp_face);
scalp_center=mean(scalp_vertex,1);
outward_test=sum(scalp_normal.*(scalp_vertex-repmat(scalp_center,size(scalp_vertex,1),1)),2);
scalp_normal(outward_test<0,:)=-scalp_normal(outward_test<0,:);

distance_to_centroid=sqrt(sum((scalp_vertex-repmat(target_centroid,size(scalp_vertex,1),1)).^2,2));
[~,seed_idx]=min(distance_to_centroid);
search_distance=sqrt(sum((scalp_vertex-repmat(scalp_vertex(seed_idx,:),size(scalp_vertex,1),1)).^2,2));
pool=find(search_distance<=search_radius_mm);
[~,order]=sort(search_distance(pool),'ascend');
pool=pool(order);
selected=zeros(0,1);
for idx=1:numel(pool)
    candidate=pool(idx);
    if(isempty(selected) || min(sqrt(sum((scalp_vertex(selected,:)-scalp_vertex(candidate,:)).^2,2)))>=grid_spacing_mm)
        selected(end+1,1)=candidate; %#ok<AGROW>
    end;
end;
if(isempty(selected))
    selected=seed_idx;
end;

selected_coord=scalp_vertex(selected,:);
selected_normal=scalp_normal(selected,:);
selected_up=zeros(size(selected_coord));
for idx=1:size(selected_coord,1)
    tangent=target_centroid-selected_coord(idx,:);
    tangent=tangent-dot(tangent,selected_normal(idx,:)).*selected_normal(idx,:);
    if(norm(tangent)<eps)
        reference=[0 0 1];
        if(abs(dot(reference,selected_normal(idx,:)))>0.9)
            reference=[0 1 0];
        end;
        tangent=reference-dot(reference,selected_normal(idx,:)).*selected_normal(idx,:);
    end;
    selected_up(idx,:)=tangent./norm(tangent);
end;

rotation_deg=double(rotation_deg(:));
if(isempty(rotation_deg))
    rotation_deg=0;
end;
n_position=size(selected_coord,1);
n_rotation=numel(rotation_deg);
n_targets=n_position*n_rotation;
target_coords=repmat(target_centroid,n_targets,1);
[target_coords_mni,native_to_mni_xfm]=etc_local_subject_surface_coords_to_mni(subjects_dir,subject,target_coords);
target_coil_rotation_deg=zeros(n_targets,1);
candidate_id=zeros(n_targets,1);
row_idx=0;
for position_idx=1:n_position
    for rotation_idx=1:n_rotation
        row_idx=row_idx+1;
        target_coil_rotation_deg(row_idx)=rotation_deg(rotation_idx);
        candidate_id(row_idx)=position_idx;
    end;
end;

scalp_center_rows=repelem(selected_coord,n_rotation,1);
scalp_normal_rows=repelem(selected_normal,n_rotation,1);
scalp_up_rows=repelem(selected_up,n_rotation,1);
target_source=repmat("tans_scalp_grid",n_targets,1);
candidate_info=struct;
candidate_info.strategy='tans_scalp_grid';
candidate_info.search_radius_mm=search_radius_mm;
candidate_info.grid_spacing_mm=grid_spacing_mm;
candidate_info.target_centroid_native=target_centroid;
candidate_info.target_centroid_mni=target_coords_mni(1,:);
candidate_info.n_roi_vertices=numel(roi_idx);
candidate_info.n_unique_coords=n_position;
candidate_info.n_rotations=n_rotation;
candidate_info.n_targets=n_targets;
candidate_info.n_atlas_targets=0;
candidate_info.n_optimized_targets=n_targets;
candidate_info.table=table(target_source,candidate_id, ...
    target_coords(:,1),target_coords(:,2),target_coords(:,3), ...
    target_coords_mni(:,1),target_coords_mni(:,2),target_coords_mni(:,3), ...
    target_coil_rotation_deg, ...
    scalp_center_rows(:,1),scalp_center_rows(:,2),scalp_center_rows(:,3), ...
    scalp_normal_rows(:,1),scalp_normal_rows(:,2),scalp_normal_rows(:,3), ...
    scalp_up_rows(:,1),scalp_up_rows(:,2),scalp_up_rows(:,3), ...
    'VariableNames',{'target_source','candidate_id', ...
    'target_native_x','target_native_y','target_native_z', ...
    'target_mni_x','target_mni_y','target_mni_z','coil_rotation_deg', ...
    'scalp_center_x','scalp_center_y','scalp_center_z', ...
    'scalp_normal_x','scalp_normal_y','scalp_normal_z', ...
    'scalp_up_x','scalp_up_y','scalp_up_z'});
fprintf('TANS scalp grid: centroid native [%1.2f %1.2f %1.2f], seed scalp vertex %d, radius=%1.1f mm, spacing=%1.1f mm, positions=%d, poses=%d.\n', ...
    target_centroid(1),target_centroid(2),target_centroid(3),seed_idx,search_radius_mm,grid_spacing_mm,n_position,n_targets);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normal=etc_local_surface_vertex_normals(vertex,face)

% read_surf returns zero-based face indices; MATLAB indexing is one-based.
face=double(face)+1;
v1=vertex(face(:,1),:);
v2=vertex(face(:,2),:);
v3=vertex(face(:,3),:);
face_normal=cross(v2-v1,v3-v1,2);
normal=zeros(size(vertex));
for face_idx=1:size(face,1)
    normal(face(face_idx,1),:)=normal(face(face_idx,1),:)+face_normal(face_idx,:);
    normal(face(face_idx,2),:)=normal(face(face_idx,2),:)+face_normal(face_idx,:);
    normal(face(face_idx,3),:)=normal(face(face_idx,3),:)+face_normal(face_idx,:);
end;
normal=normal./max(sqrt(sum(normal.^2,2)),eps);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target_coords,target_coords_mni,target_coil_rotation_deg,candidate_info,native_to_mni_xfm]=etc_local_build_target_candidates(subjects_dir, subject, label_info, decimation, max_candidate_count, min_spacing_mm, seed_strategy, rotation_deg)

roi_idx=label_info.global_vertex_index1(:);
roi_coord=double(decimation.full_vertex_coords(roi_idx,:));
roi_value=double(label_info.dfconn_value(:));
rotation_deg=double(rotation_deg(:));
if(isempty(rotation_deg))
    rotation_deg=0;
end;

switch(lower(seed_strategy))
    case 'all'
        selected_order=(1:numel(roi_idx)).';
        if(isfinite(max_candidate_count) && max_candidate_count>0 && numel(selected_order)>max_candidate_count)
            [~,sort_order]=sort(roi_value,'descend');
            selected_order=sort_order(1:max_candidate_count);
        end;
    case 'top_dfconn_farthest'
        [~,sort_order]=sort(roi_value,'descend');
        selected_order=[];
        for idx=1:numel(sort_order)
            candidate=sort_order(idx);
            if(isempty(selected_order))
                selected_order=candidate;
            else
                dist=sqrt(sum((roi_coord(selected_order,:)-roi_coord(candidate,:)).^2,2));
                if(min(dist)>=min_spacing_mm)
                    selected_order(end+1,1)=candidate; %#ok<AGROW>
                end;
            end;
            if(numel(selected_order)>=max_candidate_count)
                break;
            end;
        end;
        if(numel(selected_order)<min(max_candidate_count,numel(sort_order)))
            for idx=1:numel(sort_order)
                candidate=sort_order(idx);
                if(~ismember(candidate,selected_order))
                    selected_order(end+1,1)=candidate; %#ok<AGROW>
                end;
                if(numel(selected_order)>=max_candidate_count)
                    break;
                end;
            end;
        end;
    otherwise
        error('Unsupported candidate_seed_strategy [%s].',seed_strategy);
end;

selected_order=selected_order(:);
selected_coord=roi_coord(selected_order,:);
selected_idx=roi_idx(selected_order);
selected_value=roi_value(selected_order);
n_coord=size(selected_coord,1);
n_rot=numel(rotation_deg);
n_targets=n_coord*n_rot;
target_coords=zeros(n_targets,3);
target_coords_mni=nan(n_targets,3);
target_coil_rotation_deg=zeros(n_targets,1);
candidate_id=zeros(n_targets,1);
candidate_global_vertex_index1=zeros(n_targets,1);
candidate_dfconn_value=zeros(n_targets,1);

row_idx=0;
for coord_idx=1:n_coord
    for rot_idx=1:n_rot
        row_idx=row_idx+1;
        target_coords(row_idx,:)=selected_coord(coord_idx,:);
        target_coil_rotation_deg(row_idx)=rotation_deg(rot_idx);
        candidate_id(row_idx)=coord_idx;
        candidate_global_vertex_index1(row_idx)=selected_idx(coord_idx);
        candidate_dfconn_value(row_idx)=selected_value(coord_idx);
    end;
end;

[target_coords_mni,native_to_mni_xfm]=etc_local_subject_surface_coords_to_mni(subjects_dir, subject, target_coords);

candidate_info=struct;
candidate_info.strategy=seed_strategy;
candidate_info.max_candidate_count=max_candidate_count;
candidate_info.min_spacing_mm=min_spacing_mm;
candidate_info.n_roi_vertices=numel(roi_idx);
candidate_info.n_unique_coords=n_coord;
candidate_info.n_rotations=n_rot;
candidate_info.n_targets=n_targets;
candidate_info.n_atlas_targets=0;
candidate_info.n_optimized_targets=n_targets;
target_source=repmat("dfconn_optimized_candidate",n_targets,1);
candidate_info.table=table(target_source,candidate_id,candidate_global_vertex_index1,int32(candidate_global_vertex_index1-1), ...
    candidate_dfconn_value,target_coords_mni(:,1),target_coords_mni(:,2),target_coords_mni(:,3), ...
    target_coords(:,1),target_coords(:,2),target_coords(:,3),target_coil_rotation_deg, ...
    'VariableNames',{'target_source','candidate_id','global_vertex_index1','global_vertex_index0','dfconn_value', ...
    'target_mni_x','target_mni_y','target_mni_z','target_native_x','target_native_y','target_native_z','coil_rotation_deg'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target_coords,target_coords_mni,target_coil_rotation_deg,candidate_info]=etc_local_prepend_atlas_comparison_target(subjects_dir, subject, atlas_target_mni, atlas_rotation_deg, target_coords, target_coords_mni, target_coil_rotation_deg, candidate_info)

atlas_target_mni=double(atlas_target_mni);
if(numel(atlas_target_mni)~=3)
    error('atlas_target_coord_mni must be a 1x3 vector.');
end;
atlas_target_mni=reshape(atlas_target_mni,1,3);
atlas_native=etc_local_mni_to_subject_surface_coords(subjects_dir, subject, atlas_target_mni);
atlas_native=atlas_native(1,:);

target_coords=[atlas_native; target_coords];
target_coords_mni=[atlas_target_mni; target_coords_mni];
target_coil_rotation_deg=[double(atlas_rotation_deg); target_coil_rotation_deg(:)];

atlas_row=table("atlas_mni305",0,nan,int32(-1),nan, ...
    atlas_target_mni(1),atlas_target_mni(2),atlas_target_mni(3), ...
    atlas_native(1),atlas_native(2),atlas_native(3),double(atlas_rotation_deg), ...
    'VariableNames',candidate_info.table.Properties.VariableNames);
candidate_info.table=[atlas_row; candidate_info.table];
candidate_info.atlas_target_coord_mni=atlas_target_mni;
candidate_info.atlas_target_coord_native=atlas_native;
candidate_info.atlas_coil_rotation_deg=double(atlas_rotation_deg);
candidate_info.n_atlas_targets=1;
candidate_info.n_optimized_targets=etc_local_table_nrows(candidate_info.table)-1;
candidate_info.n_targets=etc_local_table_nrows(candidate_info.table);

fprintf('Atlas comparison target MNI305 [%1.2f %1.2f %1.2f] -> native [%1.2f %1.2f %1.2f], coil rotation %1.2f deg.\n', ...
    atlas_target_mni(1),atlas_target_mni(2),atlas_target_mni(3), ...
    atlas_native(1),atlas_native(2),atlas_native(3),double(atlas_rotation_deg));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_write_target_candidate_csv(csv_file, candidate_info)

writetable(candidate_info.table,csv_file);
fprintf('Wrote target candidate list: [%s]\n',csv_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function row=etc_local_empty_target_summary_row()

row=struct;
row.subject="";
row.target_index=nan;
row.target_source="";
row.output_stem="";
row.target_native_x=nan;
row.target_native_y=nan;
row.target_native_z=nan;
row.target_mni_x=nan;
row.target_mni_y=nan;
row.target_mni_z=nan;
row.coil_rotation_deg=nan;
row.coil_center_x=nan;
row.coil_center_y=nan;
row.coil_center_z=nan;
row.coil_center_mni_x=nan;
row.coil_center_mni_y=nan;
row.coil_center_mni_z=nan;
row.coil_axis_x=nan;
row.coil_axis_y=nan;
row.coil_axis_z=nan;
row.coil_axis_mni_x=nan;
row.coil_axis_mni_y=nan;
row.coil_axis_mni_z=nan;
row.coil_up_x=nan;
row.coil_up_y=nan;
row.coil_up_z=nan;
row.coil_up_mni_x=nan;
row.coil_up_mni_y=nan;
row.coil_up_mni_z=nan;
row.brainsight_AP_deg=nan;
row.brainsight_Lat_deg=nan;
row.brainsight_Twist_deg=nan;
row.roi_nvertices=nan;
row.roi_lh_nvertices=nan;
row.roi_rh_nvertices=nan;
row.roi_area_mm2=nan;
row.dfconn_roi_mean=nan;
row.dfconn_roi_max=nan;
row.Etotal_roi_mean_Vm=nan;
row.Etotal_roi_median_Vm=nan;
row.Etotal_roi_max_Vm=nan;
row.Enormal_abs_roi_mean_Vm=nan;
row.Etangent_roi_mean_Vm=nan;
row.TANS_OnTarget_percent=nan;
row.TANS_target_coverage_percent=nan;
row.hotspot_p99p5_on_target_percent=nan;
row.hotspot_p99p5_target_coverage_percent=nan;
row.hotspot_p99p5_overlap_area_mm2=nan;
row.hotspot_p99p5_area_mm2=nan;
row.efield_elapsed_sec=nan;
row.basis_solve_elapsed_sec=nan;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function row=etc_local_target_summary_row(subject, target_idx, output_stem_target, target_coord, target_mni, coil_rotation_deg, target_source, coil_center, coil_center_mni, coil_orientation, coil_orientation_mni, coil_up, coil_up_mni, brainsight_pose, efield, on_target, label_info)

row=etc_local_empty_target_summary_row();
row.subject=string(subject);
row.target_index=target_idx;
row.target_source=string(target_source);
row.output_stem=string(output_stem_target);
row.target_native_x=target_coord(1);
row.target_native_y=target_coord(2);
row.target_native_z=target_coord(3);
row.target_mni_x=target_mni(1);
row.target_mni_y=target_mni(2);
row.target_mni_z=target_mni(3);
row.coil_rotation_deg=coil_rotation_deg;
row.coil_center_x=coil_center(1);
row.coil_center_y=coil_center(2);
row.coil_center_z=coil_center(3);
row.coil_center_mni_x=coil_center_mni(1);
row.coil_center_mni_y=coil_center_mni(2);
row.coil_center_mni_z=coil_center_mni(3);
row.coil_axis_x=coil_orientation(1);
row.coil_axis_y=coil_orientation(2);
row.coil_axis_z=coil_orientation(3);
row.coil_axis_mni_x=coil_orientation_mni(1);
row.coil_axis_mni_y=coil_orientation_mni(2);
row.coil_axis_mni_z=coil_orientation_mni(3);
row.coil_up_x=coil_up(1);
row.coil_up_y=coil_up(2);
row.coil_up_z=coil_up(3);
row.coil_up_mni_x=coil_up_mni(1);
row.coil_up_mni_y=coil_up_mni(2);
row.coil_up_mni_z=coil_up_mni(3);
row.brainsight_AP_deg=brainsight_pose.AP_deg;
row.brainsight_Lat_deg=brainsight_pose.Lat_deg;
row.brainsight_Twist_deg=brainsight_pose.Twist_deg;
row.roi_nvertices=on_target.roi.nvertices;
row.roi_lh_nvertices=label_info.hemi(1).nvertices;
row.roi_rh_nvertices=label_info.hemi(2).nvertices;
row.roi_area_mm2=on_target.roi.area_mm2;
row.dfconn_roi_mean=label_info.dfconn_summary.mean;
row.dfconn_roi_max=label_info.dfconn_summary.max;
row.Etotal_roi_mean_Vm=on_target.Etotal.mean;
row.Etotal_roi_median_Vm=on_target.Etotal.median;
row.Etotal_roi_max_Vm=on_target.Etotal.max;
if(isfield(on_target,'Enormal_abs'))
    row.Enormal_abs_roi_mean_Vm=on_target.Enormal_abs.mean;
end;
if(isfield(on_target,'Etangent'))
    row.Etangent_roi_mean_Vm=on_target.Etangent.mean;
end;
row.TANS_OnTarget_percent=100.*on_target.TANS_OnTarget;
row.TANS_target_coverage_percent=100.*on_target.hotspot.mean_target_coverage;
row.hotspot_p99p5_on_target_percent=100.*on_target.hotspot.label.on_target;
row.hotspot_p99p5_target_coverage_percent=100.*on_target.hotspot.label.target_coverage;
if(isfield(on_target.hotspot.label,'overlap') && isfield(on_target.hotspot.label.overlap,'area_mm2'))
    row.hotspot_p99p5_overlap_area_mm2=on_target.hotspot.label.overlap.area_mm2;
elseif(isfield(on_target.hotspot.label,'overlap_area_mm2'))
    row.hotspot_p99p5_overlap_area_mm2=on_target.hotspot.label.overlap_area_mm2;
end;
row.hotspot_p99p5_area_mm2=on_target.hotspot.label.area_mm2;
if(isfield(efield,'decimation') && isfield(efield.decimation,'efield_elapsed_sec'))
    row.efield_elapsed_sec=efield.decimation.efield_elapsed_sec;
end;
if(isfield(efield,'basis_solve_elapsed_sec'))
    row.basis_solve_elapsed_sec=efield.basis_solve_elapsed_sec;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function atlas_vs_best=etc_local_make_atlas_vs_best_table(atlas_row, best_row)

metric_name=[
    "TANS_OnTarget_percent";
    "TANS_target_coverage_percent";
    "hotspot_p99p5_on_target_percent";
    "hotspot_p99p5_target_coverage_percent";
    "Etotal_roi_mean_Vm";
    "Etotal_roi_median_Vm";
    "Etotal_roi_max_Vm";
    "Enormal_abs_roi_mean_Vm";
    "Etangent_roi_mean_Vm"
    ];
n_metric=numel(metric_name);
atlas_value=nan(n_metric,1);
best_optimized_value=nan(n_metric,1);
delta_best_minus_atlas=nan(n_metric,1);
for idx=1:n_metric
    name=char(metric_name(idx));
    atlas_value(idx)=atlas_row.(name);
    best_optimized_value(idx)=best_row.(name);
    delta_best_minus_atlas(idx)=best_optimized_value(idx)-atlas_value(idx);
end;

atlas_target_native=[atlas_row.target_native_x atlas_row.target_native_y atlas_row.target_native_z];
best_target_native=[best_row.target_native_x best_row.target_native_y best_row.target_native_z];
target_distance_mm=repmat(sqrt(sum((best_target_native-atlas_target_native).^2)),n_metric,1);
atlas_rotation_deg=repmat(atlas_row.coil_rotation_deg,n_metric,1);
best_rotation_deg=repmat(best_row.coil_rotation_deg,n_metric,1);

atlas_vs_best=table(metric_name,atlas_value,best_optimized_value,delta_best_minus_atlas, ...
    target_distance_mm,atlas_rotation_deg,best_rotation_deg);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function top_result_cache=etc_local_update_top_result_cache(top_result_cache, result_entry, top_result_count)

top_result_cache{end+1}=result_entry;
score=zeros(numel(top_result_cache),4);
for idx=1:numel(top_result_cache)
    score(idx,:)=etc_local_result_score(top_result_cache{idx}.summary_row);
end;
[~,order]=sortrows(score,[-1 -2 -3 -4]);
top_result_cache=top_result_cache(order);
if(numel(top_result_cache)>top_result_count)
    top_result_cache=top_result_cache(1:top_result_count);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score=etc_local_result_score(summary_row)

score=[
    summary_row.TANS_OnTarget_percent, ...
    summary_row.TANS_target_coverage_percent, ...
    summary_row.Etotal_roi_mean_Vm, ...
    summary_row.Etotal_roi_max_Vm
    ];
score(~isfinite(score))=-inf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_save_top_result_cache(top_result_cache, root_dir, label_info, subjects_dir, flag_native_nii, flag_stc, target_coords, target_coords_mni, target_coord_mni, target_coil_rotation_deg, target_candidate_info, hotspot_label_percentile)

for idx=1:numel(top_result_cache)
    entry=top_result_cache{idx};
    efield=entry.efield;
    efield_sparse=entry.efield_sparse;
    brain_decimation=entry.brain_decimation;
    on_target=entry.on_target;
    target_info=entry.target_info;
    stc_info=entry.stc_info;
    output_stem_target=char(entry.summary_row.output_stem);

    efield_for_labels=entry.efield_ontarget;
    brain_decimation_for_labels=entry.brain_decimation_ontarget;
    if(isempty(efield_for_labels))
        efield_for_labels=efield;
        brain_decimation_for_labels=brain_decimation;
    end;
    on_target=etc_local_write_hotspot_efield_labels(on_target, efield_for_labels, brain_decimation_for_labels, root_dir, output_stem_target, hotspot_label_percentile);
    efield.on_target=on_target;
    efield_sparse.on_target=on_target;

    fn=fullfile(root_dir,sprintf('%s.mat',output_stem_target));
    fprintf('saving top e-field [%s]...\n',fn);
    if(isfield(efield,'flag_sparse') && efield.flag_sparse==0)
        save(fn, 'efield', 'efield_sparse', 'brain_decimation', 'label_info', 'on_target', 'target_info', 'target_coords', 'target_coords_mni', 'target_coord_mni', 'target_coil_rotation_deg', 'target_candidate_info', 'stc_info', 'flag_stc');
    else
        save(fn, 'efield', 'brain_decimation', 'label_info', 'on_target', 'target_info', 'target_coords', 'target_coords_mni', 'target_coord_mni', 'target_coil_rotation_deg', 'target_candidate_info', 'stc_info', 'flag_stc');
    end;

    if(flag_native_nii)
        etc_local_write_native_nii_outputs(subjects_dir, label_info.subject, efield, label_info.subject_sup_mPFC_label_files, root_dir, output_stem_target);
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_cleanup_previous_optimize_outputs(root_dir, output_stem)

patterns={
    sprintf('%s_target*',output_stem)
    sprintf('%s_all_target_metrics.csv',output_stem)
    sprintf('%s_ranked_targets.csv',output_stem)
    sprintf('%s_candidates.csv',output_stem)
    sprintf('%s_best_target.csv',output_stem)
    sprintf('%s_best_target.mat',output_stem)
    sprintf('%s_atlas_vs_best_optimized.csv',output_stem)
    };
n_deleted=0;
for pattern_idx=1:numel(patterns)
    hits=dir(fullfile(root_dir,patterns{pattern_idx}));
    for hit_idx=1:numel(hits)
        if(hits(hit_idx).isdir)
            continue;
        end;
        fn=fullfile(hits(hit_idx).folder,hits(hit_idx).name);
        delete(fn);
        n_deleted=n_deleted+1;
    end;
end;
if(n_deleted>0)
    fprintf('Deleted %d previous optimization output files for [%s].\n',n_deleted,output_stem);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf=etc_local_string_equal(value, target)

if(isstring(value))
    tf=value==string(target);
elseif(iscellstr(value))
    tf=strcmp(value,target);
elseif(ischar(value))
    tf=strcmp(cellstr(value),target);
else
    tf=strcmp(cellstr(string(value)),target);
end;
tf=tf(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=etc_local_table_nrows(T)

n=size(T,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_idx1=etc_local_read_label_local_indices(label_file, max_nvertices)

if(exist(label_file,'file')~=2)
    error('Cannot find label file [%s].',label_file);
end;
label_data=read_label('',label_file);
if(isempty(label_data))
    local_idx1=[];
else
    local_idx1=double(label_data(:,1))+1;
    local_idx1=unique(local_idx1(local_idx1>=1 & local_idx1<=max_nvertices));
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_delete_if_exists(fn)

if(exist(fn,'file')==2)
    delete(fn);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=etc_local_shell_quote(s)

s=char(s);
q=['''' strrep(s,'''','''"''"''') ''''];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_init_nav_bem_surfaces(bem_obj, head_surf_idx)

global etc_render_fsbrain

bem_idx=[1:length(bem_obj)];
bem_idx(head_surf_idx)=[];
bem_idx(end+1)=head_surf_idx;
bem_obj_tmp=bem_obj([bem_idx]);

for idx=1:length(bem_obj_tmp)
    cc=colororder;
    cc_idx=mod(idx-1,size(cc,1))+1;
    bem_obj_tmp(idx).color=cc(cc_idx,:);

    if(idx==length(bem_obj_tmp))
        bem_obj_tmp(idx).flag_show=1;
        visible='on';
    else
        bem_obj_tmp(idx).flag_show=0;
        visible='off';
    end;

    bem_obj_tmp(idx).patch_handle=patch('vertices',bem_obj_tmp(idx).vertex,'faces',bem_obj_tmp(idx).face,'edgecolor','none','facecolor',bem_obj_tmp(idx).color,'facealpha',0.2,'visible',visible);
end;

etc_render_fsbrain.surf_obj=bem_obj_tmp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_notify_nav_model_status(bem_def_status, bem_prep_status)

global etc_render_fsbrain

try
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
catch
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_apply_coil_to_nav(tms_coil_xfm_moved)

global etc_render_fsbrain

try
    if(isfield(etc_render_fsbrain,'app_tms_nav'))
        if(~isempty(etc_render_fsbrain.app_tms_nav))
            if(isvalid(etc_render_fsbrain.app_tms_nav))
                etc_tms_target_xfm_apply_nav(etc_render_fsbrain.app_tms_nav, tms_coil_xfm_moved(1:3,4).*1e3, -tms_coil_xfm_moved(1:3,3), tms_coil_xfm_moved(1:3,2));
            end;
        end;
    end;
catch
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function object_xfm=etc_local_tms_rotate_about_axis(object_xfm_base, tune_value)

% Same transform as etc_tms_target_xfm_tune(..., 4, tune_value), but without
% mutating the navigation globals/app state.
tms_coil_origin=object_xfm_base(1:3,4).*1e3;
tms_coil_axis=-object_xfm_base(1:3,3);
tms_coil_up=object_xfm_base(1:3,2);
object_xfm=object_xfm_base;

coil_center=tms_coil_origin(:);
R_c=eye(4);
R_c(1:3,4)=-coil_center(1:3)./1e3;

coil_axis=tms_coil_axis(:).';
tmp=coil_axis./norm(coil_axis);
x=tmp(1);
y=tmp(2);
z=tmp(3);

vtop=tms_coil_up(:);
vtop_perp=vtop(:)-sum(vtop(:).*tmp(:)).*tmp(:);
vup=object_xfm(1:3,2);

xx=cross(vtop_perp,vup);
c=sign(dot(xx,tms_coil_axis))*norm(xx);
deg=atan2d(c,dot(vtop_perp,vup)).*pi./180;
theta=(deg-tune_value./180*pi);

r_rotz=[cos(theta)+x*x*(1-cos(theta)), x*y*(1-cos(theta))-z*sin(theta), x*z*(1-cos(theta))+y*sin(theta);
    y*x*(1-cos(theta))+z*sin(theta) cos(theta)+y*y*(1-cos(theta)) y*z*(1-cos(theta))-x*sin(theta);
    z*x*(1-cos(theta))-y*sin(theta) z*y*(1-cos(theta))+x*sin(theta) cos(theta)+z*z*(1-cos(theta))];

R_c_inv=eye(4);
R_c_inv(1:3,4)=coil_center(1:3)./1e3;

R_rotz=eye(4);
R_rotz(1:3,1:3)=r_rotz;

object_xfm=R_c_inv*R_rotz*R_c*object_xfm;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xfm_mm=etc_local_tms_xfm_m_to_mm(xfm_m)

xfm_mm=xfm_m;
xfm_mm(1:3,4)=xfm_m(1:3,4).*1e3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target_coord_mni, target_coil_rotation_deg, target_csv_info]=etc_local_read_tms_target_csv(csv_file, subject_name, target_kind)

if(exist(csv_file,'file')~=2)
    error('Cannot find TMS target CSV [%s].',csv_file);
end;
if(nargin<3 || isempty(target_kind))
    target_kind='individualized';
end;

raw=readcell(csv_file,'Delimiter',',');
subject_name=char(string(subject_name));

switch(lower(target_kind))
    case {'atlas','generic'}
        coord_col=2:4;
        degree_col=5;
        kind='atlas';
    case {'individualized','individual','ind'}
        coord_col=8:10;
        degree_col=11;
        kind='individualized';
    otherwise
        error('Unsupported target kind [%s]. Use ''atlas'' or ''individualized''.',target_kind);
end;

target_row=[];
for row_idx=3:size(raw,1)
    row_subject=etc_local_csv_cell_to_string(raw{row_idx,1});
    if(strcmpi(row_subject, subject_name))
        target_row=row_idx;
        break;
    end;
end;

if(isempty(target_row))
    error('Subject [%s] is not listed in TMS target CSV [%s].',subject_name,csv_file);
end;

target_coord_mni=nan(1,3);
for coord_idx=1:3
    target_coord_mni(coord_idx)=etc_local_csv_cell_to_double(raw{target_row,coord_col(coord_idx)});
end;
target_coil_rotation_deg=etc_local_csv_cell_to_double(raw{target_row,degree_col});
if(isnan(target_coil_rotation_deg))
    target_coil_rotation_deg=0;
end;

if(any(isnan(target_coord_mni)))
    error('Subject [%s] has no complete %s target coordinate in [%s].',subject_name,kind,csv_file);
end;

target_csv_info=struct;
target_csv_info.file=csv_file;
target_csv_info.kind=kind;
target_csv_info.subject=subject_name;
target_csv_info.row=target_row;
target_csv_info.coord_columns=coord_col;
target_csv_info.degree_column=degree_col;
target_csv_info.coord_mni=target_coord_mni;
target_csv_info.coil_rotation_deg=target_coil_rotation_deg;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value=etc_local_csv_cell_to_string(cell_value)

if(isempty(cell_value))
    value='';
    return;
end;

try
    s=string(cell_value);
    if(any(ismissing(s)))
        value='';
    else
        value=char(strtrim(s(1)));
    end;
catch
    value='';
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value=etc_local_csv_cell_to_double(cell_value)

value=nan;
if(isempty(cell_value))
    return;
end;

if(isnumeric(cell_value) || islogical(cell_value))
    if(~isempty(cell_value))
        value=double(cell_value(1));
    end;
    return;
end;

try
    s=string(cell_value);
    if(any(ismissing(s)))
        return;
    end;
    s=strtrim(s(1));
    if(strlength(s)<1)
        return;
    end;
    value=str2double(s);
catch
    value=nan;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_stem_target=etc_local_target_output_stem(output_stem, target_idx, n_targets)

if(n_targets<=1)
    output_stem_target=output_stem;
else
    output_stem_target=sprintf('%s_target%02d',output_stem,target_idx);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function target_coords=etc_local_mni_to_subject_surface_coords(subjects_dir, subject, target_mni)

if(exist('MRIread','file')~=2)
    error('Cannot find MRIread.m on the MATLAB path; needed for MNI-to-subject target conversion.');
end;
if(exist('etc_read_xfm','file')~=2)
    error('Cannot find etc_read_xfm.m on the MATLAB path; needed to read talairach.xfm.');
end;

target_mni=double(target_mni);
if(isvector(target_mni))
    if(numel(target_mni)~=3)
        error('target_mni must be a 1x3 vector or an Nx3 matrix.');
    end;
    target_mni=reshape(target_mni,1,3);
elseif(size(target_mni,2)~=3 && size(target_mni,1)==3)
    target_mni=target_mni.';
elseif(size(target_mni,2)~=3)
    error('target_mni must be a 1x3 vector or an Nx3 matrix.');
end;

mri_file=fullfile(subjects_dir,subject,'mri','orig.mgz');
if(exist(mri_file,'file')~=2)
    error('Cannot find subject MRI [%s].',mri_file);
end;
talairach_file=fullfile(subjects_dir,subject,'mri','transforms','talairach.xfm');
if(exist(talairach_file,'file')~=2)
    error('Cannot find subject transform [%s].',talairach_file);
end;

vol=MRIread(mri_file);
vol_pre_xfm=eye(4);
if(isfield(vol,'vol_pre_xfm') && ~isempty(vol.vol_pre_xfm))
    vol_pre_xfm=vol.vol_pre_xfm;
end;
talxfm=etc_read_xfm('file_xfm',talairach_file);

target_coords=zeros(size(target_mni));
for target_idx=1:size(target_mni,1)
    target_vox=inv(vol.vox2ras)*inv(vol_pre_xfm)*inv(talxfm)*[target_mni(target_idx,:) 1].';
    target_vox=target_vox(1:3);
    surface_coord=vol.tkrvox2ras*[target_vox(:); 1];
    target_coords(target_idx,:)=surface_coord(1:3).';
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target_mni,native_to_mni_xfm]=etc_local_subject_surface_coords_to_mni(subjects_dir, subject, target_native)

if(exist('MRIread','file')~=2)
    error('Cannot find MRIread.m on the MATLAB path; needed for native-to-MNI conversion.');
end;
if(exist('etc_read_xfm','file')~=2)
    error('Cannot find etc_read_xfm.m on the MATLAB path; needed to read talairach.xfm.');
end;

target_native=double(target_native);
if(isvector(target_native))
    if(numel(target_native)~=3)
        error('target_native must be a 1x3 vector or an Nx3 matrix.');
    end;
    target_native=reshape(target_native,1,3);
elseif(size(target_native,2)~=3 && size(target_native,1)==3)
    target_native=target_native.';
elseif(size(target_native,2)~=3)
    error('target_native must be a 1x3 vector or an Nx3 matrix.');
end;

mri_file=fullfile(subjects_dir,subject,'mri','orig.mgz');
if(exist(mri_file,'file')~=2)
    error('Cannot find subject MRI [%s].',mri_file);
end;
talairach_file=fullfile(subjects_dir,subject,'mri','transforms','talairach.xfm');
if(exist(talairach_file,'file')~=2)
    error('Cannot find subject transform [%s].',talairach_file);
end;

vol=MRIread(mri_file);
vol_pre_xfm=eye(4);
if(isfield(vol,'vol_pre_xfm') && ~isempty(vol.vol_pre_xfm))
    vol_pre_xfm=vol.vol_pre_xfm;
end;
talxfm=etc_read_xfm('file_xfm',talairach_file);

% FreeSurfer surface RAS -> MNI305 transform. Keep the full affine so the
% same registration can be applied to the coil position and orientation.
native_to_mni_xfm=talxfm*vol_pre_xfm*vol.vox2ras*inv(vol.tkrvox2ras);

target_mni=zeros(size(target_native));
for target_idx=1:size(target_native,1)
    mni_coord=native_to_mni_xfm*[target_native(target_idx,:) 1].';
    target_mni(target_idx,:)=mni_coord(1:3).';
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mni_coord=etc_local_apply_native_to_mni(native_to_mni_xfm, native_coord)

native_coord=double(native_coord(:));
if(numel(native_coord)~=3)
    error('native_coord must contain exactly three coordinates.');
end;
mni_coord=(native_to_mni_xfm*[native_coord;1]).';
mni_coord=mni_coord(1:3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mni_vector=etc_local_apply_native_to_mni_vector(native_to_mni_xfm, native_vector)

native_vector=double(native_vector(:));
if(numel(native_vector)~=3)
    error('native_vector must contain exactly three components.');
end;
mni_vector=native_to_mni_xfm(1:3,1:3)*native_vector;
vector_norm=norm(mni_vector);
if(vector_norm>eps)
    mni_vector=mni_vector./vector_norm;
end;
mni_vector=mni_vector.';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pose=etc_local_compute_brainsight_pose(target_coord, coil_center, coil_up)

% Brainsight AP/Lat/Twist convention used here:
%   AP    = anterior/posterior tilt of the target-to-scalp trajectory in the
%           native sagittal (Y-Z) plane; anterior is positive.
%   Lat   = left/right tilt in the native coronal (X-Z) plane; left is
%           positive for FreeSurfer surface RAS.
%   Twist is currently reported using the native coil-up convention. Its
%   exact sign/axis mapping to the Brainsight Twist control requires
%   validation against the physical coil marker in Brainsight.
% These are the values corresponding to the AP, Lat and Twist target
% controls; the reported coil center remains the scalp position to enter as
% the crosshair origin.
target_coord=double(target_coord(:));
coil_center=double(coil_center(:));
coil_up=double(coil_up(:));

% Brainsight's AP/Lat controls describe how the trajectory leaves the target
% toward the scalp, so use target -> coil center (not coil -> target).
trajectory=coil_center-target_coord;
trajectory_norm=norm(trajectory);
if(trajectory_norm<=eps)
    error('Cannot compute Brainsight AP/Lat/Twist: target and coil center coincide.');
end;
trajectory=trajectory./trajectory_norm;

pose.AP_deg=atan2d(trajectory(2),trajectory(3));
pose.Lat_deg=atan2d(trajectory(1),trajectory(3));

% Use the trajectory rather than the scalp normal as the twist axis, since
% Brainsight twist rotates the coil around the target-to-scalp trajectory.
coil_up=coil_up-dot(coil_up,trajectory).*trajectory;
coil_up_norm=norm(coil_up);
if(coil_up_norm<=eps)
    pose.Twist_deg=nan;
    return;
end;
coil_up=coil_up./coil_up_norm;
reference_up=[0;0;1]-trajectory(3).*trajectory;
reference_up_norm=norm(reference_up);
if(reference_up_norm<=eps)
    % At a purely superior/inferior trajectory, use native anterior as the
    % zero-twist reference because projected superior is undefined.
    reference_up=[0;1;0]-trajectory(2).*trajectory;
    reference_up_norm=norm(reference_up);
end;
reference_up=reference_up./reference_up_norm;
% The simulation coil_up vector is the physical handle direction.  Negate
% the mathematical right-handed angle so that entering this value in
% Brainsight reproduces the same handle side shown by the simulation.
pose.Twist_deg=atan2d(dot(trajectory,cross(reference_up,coil_up)),dot(reference_up,coil_up));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decimation=etc_local_ensure_full_to_sparse_map(decimation, nearest_chunk_size, smooth_step)

if(nargin<3 || isempty(smooth_step))
    smooth_step=5;
end;

if(~isfield(decimation,'full_to_sparse_row_index'))
    decimation=etc_local_add_full_to_sparse_map(decimation, nearest_chunk_size, smooth_step);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_compute_on_target( ...
    efield_full, decimation, target_coord, label_files, circle_diameter_mm, circle_hemi, ...
    subjects_dir, subject, flag_exclude_sulcal_vertices, sulc_threshold, sulc_keep_direction, ...
    sulcal_empty_action, field_name, hotspot_percentile_thresholds, hotspot_hemi, ...
    flag_didt, didt_reference_A_per_us, didt_A_per_us, didt_activation_threshold_V_per_m, ...
    didt_hotspot_hemi)

% Build the target ROI and report descriptive E-field values inside it.
% The ROI can come from labels or a target-centered cortical disk; by
% default for FreeSurfer .sulc, vertices with sulc > 0 are removed so the
% retained ROI is on gyral crowns.
% This function also calls the TANS hotspot-overlap calculation below.

on_target=struct;

if(nargin<9)
    flag_exclude_sulcal_vertices=1;
end;
if(nargin<10)
    sulc_threshold=0;
end;
if(nargin<11 || isempty(sulc_keep_direction))
    sulc_keep_direction='negative';
end;
if(nargin<12 || isempty(sulcal_empty_action))
    sulcal_empty_action='nearest_gyral';
end;
if(nargin<13 || isempty(field_name))
    field_name='Etotal';
end;
if(nargin<14 || isempty(hotspot_percentile_thresholds))
    hotspot_percentile_thresholds=linspace(99.9,99,10);
end;
if(nargin<15 || isempty(hotspot_hemi))
    hotspot_hemi='both';
end;
if(nargin<16 || isempty(flag_didt))
    flag_didt=0;
end;
if(nargin<17 || isempty(didt_reference_A_per_us))
    didt_reference_A_per_us=94;
end;
if(nargin<18 || isempty(didt_A_per_us))
    didt_A_per_us=20:5:120;
end;
if(nargin<19 || isempty(didt_activation_threshold_V_per_m))
    didt_activation_threshold_V_per_m=100;
end;
if(nargin<20 || isempty(didt_hotspot_hemi))
    didt_hotspot_hemi=hotspot_hemi;
end;

if(~isempty(label_files))
    roi=etc_local_on_target_roi_from_label(decimation, label_files, subject);
else
    roi=etc_local_on_target_roi_from_circle(decimation, target_coord, circle_diameter_mm, circle_hemi);
end;

roi.vertex_index1_pre_sulcal_mask=roi.vertex_index1(:);
if(flag_exclude_sulcal_vertices)
    [sulc_value, sulc_files]=etc_local_load_sulc_values(decimation, subjects_dir, subject);
    roi=etc_local_apply_sulcal_mask_to_roi(roi, sulc_value, sulc_files, sulc_threshold, sulc_keep_direction, decimation, target_coord, sulcal_empty_action);
else
    roi.sulcal_mask=struct;
    roi.sulcal_mask.flag_exclude=0;
end;

if(isempty(roi.vertex_index1))
    if(isfield(roi,'sulcal_mask') && isfield(roi.sulcal_mask,'flag_exclude') && roi.sulcal_mask.flag_exclude)
        error('On-target ROI is empty after sulcal exclusion. Initial ROI had %d vertices; none survived %s %1.3f. For circular ROIs, set on_target_empty_after_sulcal_mask_action=''nearest_gyral'' to recenter, increase on_target_circle_diameter_mm, or set on_target_exclude_sulcal_vertices=0.', ...
            roi.sulcal_mask.nvertices_before, roi.sulcal_mask.rule, roi.sulcal_mask.threshold);
    else
        error('On-target ROI is empty.');
    end;
end;

vertex_area=etc_local_vertex_surface_area(decimation);

roi.nvertices=length(roi.vertex_index1);
roi.vertex_index0=int32(roi.vertex_index1(:)-1);
roi.center_target_coord=target_coord(:).';
roi.vertex_coords=decimation.full_vertex_coords(roi.vertex_index1,:);
roi.vertex_area_mm2=vertex_area(roi.vertex_index1);
roi.area_mm2=sum(roi.vertex_area_mm2);
roi.hemi_count=etc_local_count_vertices_by_hemi(decimation, roi.vertex_index1);

on_target.roi=roi;

on_target.Etotal=etc_local_value_summary(efield_full.Etotal(roi.vertex_index1));
if(isfield(efield_full,'Enormal'))
    on_target.Enormal=etc_local_value_summary(efield_full.Enormal(roi.vertex_index1));
    on_target.Enormal_abs=etc_local_value_summary(abs(efield_full.Enormal(roi.vertex_index1)));
end;
if(isfield(efield_full,'Etangent'))
    on_target.Etangent=etc_local_value_summary(efield_full.Etangent(roi.vertex_index1));
end;
if(isfield(efield_full,'E'))
    E_roi=efield_full.E(roi.vertex_index1,:);
    on_target.E.mean_vector=mean(E_roi,1);
    on_target.E.mean_vector_norm=norm(on_target.E.mean_vector);
end;

[~,center_rel_idx]=min(sum((roi.vertex_coords-repmat(target_coord(:).',size(roi.vertex_coords,1),1)).^2,2));
center_vertex_idx=roi.vertex_index1(center_rel_idx);
on_target.target_vertex.index1=center_vertex_idx;
on_target.target_vertex.index0=int32(center_vertex_idx-1);
on_target.target_vertex.coord=decimation.full_vertex_coords(center_vertex_idx,:);
on_target.target_vertex.distance_to_target_mm=sqrt(sum((on_target.target_vertex.coord-target_coord(:).').^2));
on_target.target_vertex.Etotal=efield_full.Etotal(center_vertex_idx);
if(isfield(efield_full,'Enormal'))
    on_target.target_vertex.Enormal=efield_full.Enormal(center_vertex_idx);
end;
if(isfield(efield_full,'Etangent'))
    on_target.target_vertex.Etangent=efield_full.Etangent(center_vertex_idx);
end;
if(isfield(efield_full,'E'))
    on_target.target_vertex.E=efield_full.E(center_vertex_idx,:);
end;

% Keep this as a conventional dose/strength summary inside the ROI.
% It is not the TANS on-target specificity metric.
on_target.E_ROI=on_target.Etotal.mean;

% The TANS metric is the surface-area weighted fraction of the E-field
% magnitude hotspot that lies inside the ROI, averaged across the configured
% high-percentile hotspot thresholds.
on_target=etc_local_compute_tans_hotspot_score(on_target, efield_full, decimation, vertex_area, field_name, hotspot_percentile_thresholds, hotspot_hemi);
if(flag_didt)
    on_target=etc_local_compute_didt_dose_on_target( ...
        on_target, efield_full, decimation, vertex_area, field_name, ...
        didt_reference_A_per_us, didt_A_per_us, ...
        didt_activation_threshold_V_per_m, didt_hotspot_hemi);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sulc_value, sulc_files]=etc_local_load_sulc_values(decimation, subjects_dir, subject)

sulc_value=nan(decimation.full_nvertices,1);
sulc_files=cell(length(decimation.hemi),1);

for hemi_idx=1:length(decimation.hemi)
    hemi_name=lower(decimation.hemi(hemi_idx).name);
    sulc_file=fullfile(subjects_dir,subject,'surf',sprintf('%s.sulc',hemi_name));
    if(exist(sulc_file,'file')~=2)
        error('Cannot find sulcal-depth file [%s].',sulc_file);
    end;
    if(exist('read_curv','file')~=2)
        error('Cannot find read_curv.m on the MATLAB path; needed to read [%s].',sulc_file);
    end;

    hemi_sulc=read_curv(sulc_file);
    hemi_sulc=hemi_sulc(:);
    full_count=decimation.hemi(hemi_idx).full_vertex_count;
    if(length(hemi_sulc)~=full_count)
        error('Sulc file [%s] has %d vertices, expected %d.',sulc_file,length(hemi_sulc),full_count);
    end;

    full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
    full_idx=(full_offset+1):(full_offset+full_count);
    sulc_value(full_idx)=hemi_sulc;
    sulc_files{hemi_idx}=sulc_file;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi=etc_local_apply_sulcal_mask_to_roi(roi, sulc_value, sulc_files, sulc_threshold, sulc_keep_direction, decimation, target_coord, empty_action)

if(nargin<5 || isempty(sulc_keep_direction))
    sulc_keep_direction='negative';
end;
if(nargin<8 || isempty(empty_action))
    empty_action='nearest_gyral';
end;

roi_idx=roi.vertex_index1(:);
roi_sulc=sulc_value(roi_idx);
keep_idx=etc_local_sulcal_keep_index(roi_sulc, sulc_threshold, sulc_keep_direction);

roi.vertex_index1=roi_idx(keep_idx);
roi.sulcal_mask=struct;
roi.sulcal_mask.flag_exclude=1;
roi.sulcal_mask.rule=etc_local_sulcal_rule_string(sulc_keep_direction);
roi.sulcal_mask.keep_direction=sulc_keep_direction;
roi.sulcal_mask.threshold=sulc_threshold;
roi.sulcal_mask.files=sulc_files;
roi.sulcal_mask.nvertices_before=length(roi_idx);
roi.sulcal_mask.nvertices_after=length(roi.vertex_index1);
roi.sulcal_mask.nvertices_after_pre_recenter=length(roi.vertex_index1);
roi.sulcal_mask.nvertices_removed=length(roi_idx)-length(roi.vertex_index1);
roi.sulcal_mask.percent_after=100.*roi.sulcal_mask.nvertices_after./max(1,roi.sulcal_mask.nvertices_before);
roi.sulcal_mask.percent_after_pre_recenter=100.*roi.sulcal_mask.nvertices_after_pre_recenter./max(1,roi.sulcal_mask.nvertices_before);
roi.sulcal_mask.removed_vertex_index1=roi_idx(setdiff(1:length(roi_idx),keep_idx));
roi.sulcal_mask.empty_action=empty_action;
roi.sulcal_mask.flag_recentered=0;

if(isempty(roi.vertex_index1))
    if(strcmpi(empty_action,'nearest_gyral') && strcmpi(roi.method,'circle'))
        roi=etc_local_recenter_empty_circle_roi_to_gyral(roi, sulc_value, sulc_threshold, sulc_keep_direction, decimation, target_coord);
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi=etc_local_recenter_empty_circle_roi_to_gyral(roi, sulc_value, sulc_threshold, sulc_keep_direction, decimation, target_coord)

candidate_idx=etc_local_full_indices_for_hemi(decimation, roi.circle_hemi, target_coord);
candidate_sulc=sulc_value(candidate_idx);
gyral_idx=candidate_idx(etc_local_sulcal_keep_index(candidate_sulc, sulc_threshold, sulc_keep_direction));

if(isempty(gyral_idx))
    return;
end;

gyral_coord=decimation.full_vertex_coords(gyral_idx,:);
dist_to_target=sqrt(sum((gyral_coord-repmat(target_coord(:).',size(gyral_coord,1),1)).^2,2));
[min_dist,center_rel_idx]=min(dist_to_target);
center_vertex_idx=gyral_idx(center_rel_idx);
center_coord=decimation.full_vertex_coords(center_vertex_idx,:);

candidate_coord=decimation.full_vertex_coords(candidate_idx,:);
dist_to_center=sqrt(sum((candidate_coord-repmat(center_coord,size(candidate_coord,1),1)).^2,2));
candidate_sulc=sulc_value(candidate_idx);
keep_idx=etc_local_sulcal_keep_index(candidate_sulc, sulc_threshold, sulc_keep_direction);
roi_idx=candidate_idx(intersect(find(dist_to_center<=roi.circle_radius_mm),keep_idx));

if(isempty(roi_idx))
    roi_idx=center_vertex_idx;
end;

roi.vertex_index1=roi_idx(:);
roi.distance_metric='euclidean_from_nearest_gyral_vertex_to_surface_vertices_after_sulcal_mask';
roi.center_vertex_index1=center_vertex_idx;
roi.center_vertex_index0=int32(center_vertex_idx-1);
roi.center_vertex_coord=center_coord;
roi.center_distance_to_target_mm=min_dist;
roi.sulcal_mask.flag_recentered=1;
roi.sulcal_mask.recenter_vertex_index1=center_vertex_idx;
roi.sulcal_mask.recenter_vertex_index0=int32(center_vertex_idx-1);
roi.sulcal_mask.recenter_vertex_coord=center_coord;
roi.sulcal_mask.recenter_distance_to_target_mm=min_dist;
roi.sulcal_mask.nvertices_after=length(roi.vertex_index1);
roi.sulcal_mask.percent_after=100.*roi.sulcal_mask.nvertices_after./max(1,roi.sulcal_mask.nvertices_before);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keep_idx=etc_local_sulcal_keep_index(sulc_value, sulc_threshold, sulc_keep_direction)

switch(lower(sulc_keep_direction))
    case {'negative','neg','le','less_equal'}
        keep_idx=find(isfinite(sulc_value) & sulc_value<=sulc_threshold);
    case {'positive','pos','ge','greater_equal'}
        keep_idx=find(isfinite(sulc_value) & sulc_value>=sulc_threshold);
    otherwise
        error('Unsupported sulc_keep_direction [%s]. Use ''negative'' or ''positive''.',sulc_keep_direction);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rule=etc_local_sulcal_rule_string(sulc_keep_direction)

switch(lower(sulc_keep_direction))
    case {'negative','neg','le','less_equal'}
        rule='sulc_value <= threshold';
    case {'positive','pos','ge','greater_equal'}
        rule='sulc_value >= threshold';
    otherwise
        rule=sprintf('unsupported sulc_keep_direction [%s]',sulc_keep_direction);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vertex_area=etc_local_vertex_surface_area(decimation)

vertex_area=zeros(decimation.full_nvertices,1);
faces=double(decimation.full_faces);
if(isempty(faces))
    return;
end;
if(min(faces(:))==0)
    faces=faces+1;
end;

valid_face=find(all(faces>=1 & faces<=decimation.full_nvertices,2));
faces=faces(valid_face,:);
if(isempty(faces))
    return;
end;

vertex=double(decimation.full_vertex_coords);
v1=vertex(faces(:,1),:);
v2=vertex(faces(:,2),:);
v3=vertex(faces(:,3),:);
face_area=0.5.*sqrt(sum(cross(v2-v1,v3-v1,2).^2,2));

vertex_area=accumarray(faces(:),repmat(face_area./3,3,1),[decimation.full_nvertices 1],@sum,0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_compute_tans_hotspot_score(on_target, efield_full, decimation, vertex_area, field_name, percentile_thresholds, hotspot_hemi)

% Compute the paper-style TANS on-target score. The TANS code uses SimNIBS
% magnE; in this script the matching field is Etotal. For each percentile
% threshold, the hotspot is Etotal > prctile(Etotal,pct), and on-target is:
%   area(hotspot intersect ROI) / area(hotspot)
% target_coverage is also stored, but it is the converse diagnostic:
%   area(hotspot intersect ROI) / area(ROI)

if(~isfield(efield_full,field_name))
    error('E-field structure does not contain field [%s].',field_name);
end;

Efield=efield_full.(field_name);
if(size(Efield,2)>1)
    Efield=sqrt(sum(Efield.^2,2));
end;
Efield=Efield(:);

candidate_idx=etc_local_full_indices_for_hemi(decimation, hotspot_hemi, on_target.roi.center_target_coord);
candidate_value=Efield(candidate_idx);
finite_idx=find(isfinite(candidate_value));
if(isempty(finite_idx))
    error('No finite %s values available for TANS hotspot calculation.',field_name);
end;

target_mask=false(decimation.full_nvertices,1);
target_mask(on_target.roi.vertex_index1)=true;
target_area=sum(vertex_area(on_target.roi.vertex_index1));

percentile_thresholds=percentile_thresholds(:).';
threshold=struct([]);
score=zeros(length(percentile_thresholds),1);
target_coverage=zeros(length(percentile_thresholds),1);
threshold_value=zeros(length(percentile_thresholds),1);

finite_candidate_value=candidate_value(finite_idx);
for pct_idx=1:length(percentile_thresholds)
    pct=percentile_thresholds(pct_idx);
    threshold_value(pct_idx)=etc_local_percentile(finite_candidate_value,pct);

    hotspot_idx=candidate_idx(find(candidate_value>threshold_value(pct_idx) & isfinite(candidate_value)));
    overlap_idx=hotspot_idx(find(target_mask(hotspot_idx)));

    hotspot_area=sum(vertex_area(hotspot_idx));
    overlap_area=sum(vertex_area(overlap_idx));

    if(hotspot_area>0)
        score(pct_idx)=overlap_area./hotspot_area;
    else
        score(pct_idx)=0;
    end;
    if(target_area>0)
        target_coverage(pct_idx)=overlap_area./target_area;
    else
        target_coverage(pct_idx)=0;
    end;

    threshold(pct_idx).percentile=pct;
    threshold(pct_idx).threshold_value=threshold_value(pct_idx);
    threshold(pct_idx).vertex_index1=hotspot_idx(:);
    threshold(pct_idx).vertex_index0=int32(hotspot_idx(:)-1);
    threshold(pct_idx).nvertices=length(hotspot_idx);
    threshold(pct_idx).area_mm2=hotspot_area;
    threshold(pct_idx).overlap_vertex_index1=overlap_idx(:);
    threshold(pct_idx).overlap_vertex_index0=int32(overlap_idx(:)-1);
    threshold(pct_idx).overlap_nvertices=length(overlap_idx);
    threshold(pct_idx).overlap_area_mm2=overlap_area;
    threshold(pct_idx).on_target=score(pct_idx);
    threshold(pct_idx).target_coverage=target_coverage(pct_idx);
end;

hotspot=struct;
hotspot.method='TANS_surface_area_weighted_hotspot_overlap';
hotspot.field_name=field_name;
hotspot.percentile_thresholds=percentile_thresholds;
hotspot.threshold_values=threshold_value(:).';
hotspot.threshold=threshold;
hotspot.hotspot_hemi=hotspot_hemi;
hotspot.target_area_mm2=target_area;
hotspot.on_target_by_threshold=score(:).';
hotspot.target_coverage_by_threshold=target_coverage(:).';
hotspot.mean_on_target=mean(score);
hotspot.mean_target_coverage=mean(target_coverage);

on_target.hotspot=hotspot;
on_target.TANS_OnTarget=hotspot.mean_on_target;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_compute_didt_dose_on_target( ...
    on_target, efield_full, decimation, vertex_area, field_name, ...
    didt_reference_A_per_us, didt_A_per_us, activation_threshold_V_per_m, hotspot_hemi)

% For a fixed coil pose, the E-field is linear in dI/dt. This computes the
% absolute-threshold TANS dose curve without rerunning the BEM solve:
%   E(dI/dt) = E(reference) * dI/dt / reference_dI/dt
% For each activation threshold, on-target is:
%   area(suprathreshold hotspot intersect ROI) / area(suprathreshold hotspot)
% target_coverage is the converse diagnostic:
%   area(suprathreshold hotspot intersect ROI) / area(ROI)

if(~isfield(efield_full,field_name))
    error('E-field structure does not contain field [%s].',field_name);
end;
if(~isfinite(didt_reference_A_per_us) || didt_reference_A_per_us<=0)
    error('didt_reference_A_per_us must be positive and finite.');
end;

didt_A_per_us=double(didt_A_per_us(:)).';
activation_threshold_V_per_m=double(activation_threshold_V_per_m(:)).';
if(isempty(didt_A_per_us) || any(~isfinite(didt_A_per_us)) || any(didt_A_per_us<=0))
    error('didt_A_per_us must contain positive finite values.');
end;
if(isempty(activation_threshold_V_per_m) || any(~isfinite(activation_threshold_V_per_m)) || any(activation_threshold_V_per_m<=0))
    error('activation_threshold_V_per_m must contain positive finite values.');
end;

Efield=efield_full.(field_name);
if(size(Efield,2)>1)
    Efield=sqrt(sum(Efield.^2,2));
end;
Efield=Efield(:);

candidate_idx=etc_local_full_indices_for_hemi(decimation, hotspot_hemi, on_target.roi.center_target_coord);
candidate_value_ref=Efield(candidate_idx);
finite_idx=find(isfinite(candidate_value_ref));
if(isempty(finite_idx))
    error('No finite %s values available for dI/dt dose calculation.',field_name);
end;
candidate_idx=candidate_idx(finite_idx);
candidate_value_ref=candidate_value_ref(finite_idx);

target_mask=false(decimation.full_nvertices,1);
target_mask(on_target.roi.vertex_index1)=true;
target_area=sum(vertex_area(on_target.roi.vertex_index1));

n_threshold=length(activation_threshold_V_per_m);
n_didt=length(didt_A_per_us);
scale=didt_A_per_us./didt_reference_A_per_us;

on_target_fraction=zeros(n_threshold,n_didt);
target_coverage_fraction=zeros(n_threshold,n_didt);
hotspot_area_mm2=zeros(n_threshold,n_didt);
overlap_area_mm2=zeros(n_threshold,n_didt);
hotspot_nvertices=zeros(n_threshold,n_didt);
overlap_nvertices=zeros(n_threshold,n_didt);
effective_reference_threshold_V_per_m=zeros(n_threshold,n_didt);
threshold=struct([]);

for threshold_idx=1:n_threshold
    threshold_value=activation_threshold_V_per_m(threshold_idx);
    entry=struct([]);

    for didt_idx=1:n_didt
        effective_reference_threshold_V_per_m(threshold_idx,didt_idx)=threshold_value./scale(didt_idx);
        hotspot_idx=candidate_idx(find(candidate_value_ref.*scale(didt_idx)>=threshold_value));
        overlap_idx=hotspot_idx(find(target_mask(hotspot_idx)));

        hotspot_area=sum(vertex_area(hotspot_idx));
        overlap_area=sum(vertex_area(overlap_idx));

        hotspot_area_mm2(threshold_idx,didt_idx)=hotspot_area;
        overlap_area_mm2(threshold_idx,didt_idx)=overlap_area;
        hotspot_nvertices(threshold_idx,didt_idx)=length(hotspot_idx);
        overlap_nvertices(threshold_idx,didt_idx)=length(overlap_idx);

        if(hotspot_area>0)
            on_target_fraction(threshold_idx,didt_idx)=overlap_area./hotspot_area;
        else
            on_target_fraction(threshold_idx,didt_idx)=0;
        end;
        if(target_area>0)
            target_coverage_fraction(threshold_idx,didt_idx)=overlap_area./target_area;
        else
            target_coverage_fraction(threshold_idx,didt_idx)=0;
        end;

        entry(didt_idx).didt_A_per_us=didt_A_per_us(didt_idx);
        entry(didt_idx).scale=scale(didt_idx);
        entry(didt_idx).activation_threshold_V_per_m=threshold_value;
        entry(didt_idx).effective_reference_threshold_V_per_m=effective_reference_threshold_V_per_m(threshold_idx,didt_idx);
        entry(didt_idx).hotspot_nvertices=hotspot_nvertices(threshold_idx,didt_idx);
        entry(didt_idx).overlap_nvertices=overlap_nvertices(threshold_idx,didt_idx);
        entry(didt_idx).hotspot_area_mm2=hotspot_area;
        entry(didt_idx).overlap_area_mm2=overlap_area;
        entry(didt_idx).on_target=on_target_fraction(threshold_idx,didt_idx);
        entry(didt_idx).target_coverage=target_coverage_fraction(threshold_idx,didt_idx);
    end;

    valid_hotspot=find(hotspot_area_mm2(threshold_idx,:)>0);
    if(isempty(valid_hotspot))
        best_idx=1;
    else
        [~,best_rel_idx]=max(on_target_fraction(threshold_idx,valid_hotspot));
        best_idx=valid_hotspot(best_rel_idx);
    end;

    threshold(threshold_idx).activation_threshold_V_per_m=threshold_value;
    threshold(threshold_idx).entry=entry;
    threshold(threshold_idx).on_target_by_didt=on_target_fraction(threshold_idx,:);
    threshold(threshold_idx).target_coverage_by_didt=target_coverage_fraction(threshold_idx,:);
    threshold(threshold_idx).hotspot_area_mm2=hotspot_area_mm2(threshold_idx,:);
    threshold(threshold_idx).overlap_area_mm2=overlap_area_mm2(threshold_idx,:);
    threshold(threshold_idx).hotspot_nvertices=hotspot_nvertices(threshold_idx,:);
    threshold(threshold_idx).overlap_nvertices=overlap_nvertices(threshold_idx,:);
    threshold(threshold_idx).effective_reference_threshold_V_per_m=effective_reference_threshold_V_per_m(threshold_idx,:);
    threshold(threshold_idx).best_index=best_idx;
    threshold(threshold_idx).best_didt_A_per_us=didt_A_per_us(best_idx);
    threshold(threshold_idx).best_on_target=on_target_fraction(threshold_idx,best_idx);
    threshold(threshold_idx).best_target_coverage=target_coverage_fraction(threshold_idx,best_idx);
    threshold(threshold_idx).best_hotspot_area_mm2=hotspot_area_mm2(threshold_idx,best_idx);
    threshold(threshold_idx).best_overlap_area_mm2=overlap_area_mm2(threshold_idx,best_idx);
end;

roi_value_ref=Efield(on_target.roi.vertex_index1);
roi_value_ref=roi_value_ref(isfinite(roi_value_ref));
if(isempty(roi_value_ref))
    roi_mean_ref=nan;
    roi_max_ref=nan;
else
    roi_mean_ref=mean(roi_value_ref);
    roi_max_ref=max(roi_value_ref);
end;

dose=struct;
dose.method='absolute_threshold_didt_scaled_hotspot_overlap';
dose.field_name=field_name;
dose.hotspot_hemi=hotspot_hemi;
dose.reference_A_per_us=didt_reference_A_per_us;
dose.didt_A_per_us=didt_A_per_us;
dose.scale=scale;
dose.activation_threshold_V_per_m=activation_threshold_V_per_m;
dose.target_area_mm2=target_area;
dose.reference_max_E_V_per_m=max(candidate_value_ref);
dose.max_E_by_didt_V_per_m=max(candidate_value_ref).*scale;
dose.roi_mean_E_by_didt_V_per_m=roi_mean_ref.*scale;
dose.roi_max_E_by_didt_V_per_m=roi_max_ref.*scale;
dose.effective_reference_threshold_V_per_m=effective_reference_threshold_V_per_m;
dose.on_target_fraction=on_target_fraction;
dose.target_coverage_fraction=target_coverage_fraction;
dose.hotspot_area_mm2=hotspot_area_mm2;
dose.overlap_area_mm2=overlap_area_mm2;
dose.hotspot_nvertices=hotspot_nvertices;
dose.overlap_nvertices=overlap_nvertices;
dose.threshold=threshold;

on_target.didt=dose;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_write_didt_dose_csv(on_target, output_dir, output_stem)

if(~isfield(on_target,'didt'))
    return;
end;
if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end;

dose=on_target.didt;
fn=fullfile(output_dir,sprintf('%s_didt_on_target.csv',output_stem));
fid=fopen(fn,'w');
if(fid<0)
    error('Cannot open [%s] for writing.',fn);
end;

target_coord=on_target.roi.center_target_coord;
fprintf(fid,[ ...
    'target_x_mm,target_y_mm,target_z_mm,field_name,reference_didt_A_per_us,' ...
    'didt_A_per_us,scale,activation_threshold_V_per_m,effective_reference_threshold_V_per_m,' ...
    'on_target_percent,target_coverage_percent,hotspot_area_mm2,overlap_area_mm2,' ...
    'hotspot_nvertices,overlap_nvertices,max_E_at_didt_V_per_m,' ...
    'roi_mean_E_at_didt_V_per_m,roi_max_E_at_didt_V_per_m\n']);
for threshold_idx=1:length(dose.activation_threshold_V_per_m)
    for didt_idx=1:length(dose.didt_A_per_us)
        fprintf(fid,'%1.6f,%1.6f,%1.6f,%s,%1.6f,%1.6f,%1.9f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%d,%d,%1.6f,%1.6f,%1.6f\n', ...
            target_coord(1), target_coord(2), target_coord(3), dose.field_name, dose.reference_A_per_us, ...
            dose.didt_A_per_us(didt_idx), dose.scale(didt_idx), dose.activation_threshold_V_per_m(threshold_idx), ...
            dose.effective_reference_threshold_V_per_m(threshold_idx,didt_idx), ...
            100.*dose.on_target_fraction(threshold_idx,didt_idx), ...
            100.*dose.target_coverage_fraction(threshold_idx,didt_idx), ...
            dose.hotspot_area_mm2(threshold_idx,didt_idx), dose.overlap_area_mm2(threshold_idx,didt_idx), ...
            dose.hotspot_nvertices(threshold_idx,didt_idx), dose.overlap_nvertices(threshold_idx,didt_idx), ...
            dose.max_E_by_didt_V_per_m(didt_idx), dose.roi_mean_E_by_didt_V_per_m(didt_idx), dose.roi_max_E_by_didt_V_per_m(didt_idx));
    end;
end;
fclose(fid);

on_target.didt.csv_file=fn;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_print_didt_dose_on_target(on_target)

if(~isfield(on_target,'didt'))
    return;
end;

dose=on_target.didt;
fprintf('dI/dt dose on-target using %s, reference=%1.3g A/us, hotspot hemi=%s\n', ...
    dose.field_name, dose.reference_A_per_us, dose.hotspot_hemi);
if(isfield(dose,'csv_file'))
    fprintf('dI/dt dose CSV: [%s]\n',dose.csv_file);
end;

for threshold_idx=1:length(dose.threshold)
    threshold=dose.threshold(threshold_idx);
    best_idx=threshold.best_index;
    fprintf('  activation threshold %1.3g V/m: best on-target=%1.2f%% at dI/dt=%1.3g A/us; target coverage=%1.2f%%; hotspot area=%1.1f mm^2\n', ...
        threshold.activation_threshold_V_per_m, 100.*threshold.best_on_target, threshold.best_didt_A_per_us, ...
        100.*threshold.best_target_coverage, threshold.best_hotspot_area_mm2);
    fprintf('    dI/dt A/us: ');
    fprintf('%6.1f',dose.didt_A_per_us);
    fprintf('\n');
    fprintf('    on-target%%: ');
    fprintf('%6.1f',100.*dose.on_target_fraction(threshold_idx,:));
    fprintf('\n');
    fprintf('    targetcov%%: ');
    fprintf('%6.1f',100.*dose.target_coverage_fraction(threshold_idx,:));
    fprintf('\n');
    fprintf('    hotspot n : ');
    fprintf('%6d',dose.hotspot_nvertices(threshold_idx,:));
    fprintf('\n');
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi=etc_local_on_target_roi_from_circle(decimation, target_coord, circle_diameter_mm, circle_hemi)

radius_mm=circle_diameter_mm./2;
candidate_idx=etc_local_full_indices_for_hemi(decimation, circle_hemi, target_coord);

candidate_coord=decimation.full_vertex_coords(candidate_idx,:);
[~,center_rel_idx]=min(sum((candidate_coord-repmat(target_coord(:).',size(candidate_coord,1),1)).^2,2));
center_vertex_idx=candidate_idx(center_rel_idx);
center_coord=decimation.full_vertex_coords(center_vertex_idx,:);

dist_mm=sqrt(sum((candidate_coord-repmat(target_coord(:).',size(candidate_coord,1),1)).^2,2));
roi_idx=candidate_idx(find(dist_mm<=radius_mm));

if(isempty(roi_idx))
    roi_idx=center_vertex_idx;
end;

roi=struct;
roi.method='circle';
roi.circle_diameter_mm=circle_diameter_mm;
roi.circle_radius_mm=radius_mm;
roi.circle_hemi=circle_hemi;
roi.distance_metric='euclidean_from_target_coord_to_surface_vertices';
roi.center_vertex_index1=center_vertex_idx;
roi.center_vertex_index0=int32(center_vertex_idx-1);
roi.center_vertex_coord=center_coord;
roi.center_distance_to_target_mm=sqrt(sum((center_coord-target_coord(:).').^2));
roi.vertex_index1=roi_idx(:);
roi.hemi_count=etc_local_count_vertices_by_hemi(decimation, roi.vertex_index1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi=etc_local_on_target_roi_from_label(decimation, label_files, subject)

if(ischar(label_files))
    label_files={label_files};
end;

roi_idx=[];
roi_label=struct([]);

for label_idx=1:length(label_files)
    label_file=label_files{label_idx};
    hemi=etc_local_infer_hemi_from_label_file(label_file);
    if(isempty(hemi))
        error('Cannot infer hemisphere from label file [%s]. Use names containing lh or rh.',label_file);
    end;

    label_data=etc_local_read_label_file(subject, label_file);
    hemi_idx=find(strcmpi({decimation.hemi.name},hemi),1);
    if(isempty(hemi_idx))
        error('Cannot find hemisphere [%s] in decimation metadata.',hemi);
    end;

    local_idx1=double(label_data(:,1))+1;
    local_idx1=local_idx1(local_idx1>=1 & local_idx1<=decimation.hemi(hemi_idx).full_vertex_count);
    global_idx1=local_idx1+decimation.hemi(hemi_idx).full_vertex_offset;
    roi_idx=cat(1,roi_idx(:),global_idx1(:));

    roi_label(label_idx).file=label_file;
    roi_label(label_idx).hemi=hemi;
    roi_label(label_idx).nvertices=length(global_idx1);
end;

roi=struct;
roi.method='label';
roi.label=roi_label;
roi.vertex_index1=unique(roi_idx(:));
roi.hemi_count=etc_local_count_vertices_by_hemi(decimation, roi.vertex_index1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_write_on_target_circle_labels(on_target, decimation, output_dir, output_stem)

if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end;

diameter_str=sprintf('%g',on_target.roi.circle_diameter_mm);
diameter_str=strrep(diameter_str,'.','p');
label_prefix=fullfile(output_dir,sprintf('%s_circle%smm',output_stem,diameter_str));
label_files=etc_local_write_vertex_labels_by_hemi(decimation, on_target.roi.vertex_index1, ones(size(on_target.roi.vertex_index1)), label_prefix);

on_target.roi.circle_label_files=label_files;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_write_hotspot_efield_labels(on_target, efield_full, decimation, output_dir, output_stem, label_percentile)

if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end;

if(~isfield(on_target,'hotspot'))
    error('TANS hotspot score must be computed before writing hotspot labels.');
end;
if(isempty(label_percentile))
    label_percentile=99.5;
end;

percentiles=[on_target.hotspot.threshold.percentile];
[~,threshold_idx]=min(abs(percentiles-label_percentile));
threshold=on_target.hotspot.threshold(threshold_idx);

Efield=efield_full.(on_target.hotspot.field_name);
if(size(Efield,2)>1)
    Efield=sqrt(sum(Efield.^2,2));
end;
Efield=Efield(:);

threshold_tag=sprintf('p%g',threshold.percentile);
threshold_tag=strrep(threshold_tag,'.','p');

hotspot_idx=threshold.vertex_index1(:);
hotspot_value=Efield(hotspot_idx);
label_prefix=fullfile(output_dir,sprintf('%s_hotspot_%s_%s',output_stem,on_target.hotspot.field_name,threshold_tag));
label_files=etc_local_write_vertex_labels_by_hemi(decimation, hotspot_idx, hotspot_value, label_prefix);

overlap_idx=threshold.overlap_vertex_index1(:);
overlap_value=Efield(overlap_idx);
overlap_label_prefix=fullfile(output_dir,sprintf('%s_hotspot_%s_%s_overlap_target',output_stem,on_target.hotspot.field_name,threshold_tag));
overlap_label_files=etc_local_write_vertex_labels_by_hemi(decimation, overlap_idx, overlap_value, overlap_label_prefix);

label=struct;
label.percentile=threshold.percentile;
label.threshold_value=threshold.threshold_value;
label.vertex_index1=hotspot_idx;
label.vertex_index0=int32(hotspot_idx-1);
label.nvertices=length(hotspot_idx);
label.area_mm2=threshold.area_mm2;
label.on_target=threshold.on_target;
label.target_coverage=threshold.target_coverage;
label.files=label_files;
label.overlap.vertex_index1=overlap_idx;
label.overlap.vertex_index0=int32(overlap_idx-1);
label.overlap.nvertices=length(overlap_idx);
label.overlap.area_mm2=threshold.overlap_area_mm2;
label.overlap.files=overlap_label_files;

on_target.hotspot.label_percentile=threshold.percentile;
on_target.hotspot.label_threshold_index=threshold_idx;
on_target.hotspot.label=label;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function on_target=etc_local_attach_hotspot_label_summary(on_target, label_percentile)

if(~isfield(on_target,'hotspot'))
    error('TANS hotspot score must be computed before attaching hotspot summary.');
end;
if(isempty(label_percentile))
    label_percentile=99.5;
end;

percentiles=[on_target.hotspot.threshold.percentile];
[~,threshold_idx]=min(abs(percentiles-label_percentile));
threshold=on_target.hotspot.threshold(threshold_idx);

label=struct;
label.percentile=threshold.percentile;
label.threshold_value=threshold.threshold_value;
label.vertex_index1=threshold.vertex_index1(:);
label.vertex_index0=threshold.vertex_index0(:);
label.nvertices=threshold.nvertices;
label.area_mm2=threshold.area_mm2;
label.on_target=threshold.on_target;
label.target_coverage=threshold.target_coverage;
label.files={};
label.overlap.vertex_index1=threshold.overlap_vertex_index1(:);
label.overlap.vertex_index0=threshold.overlap_vertex_index0(:);
label.overlap.nvertices=threshold.overlap_nvertices;
label.overlap.area_mm2=threshold.overlap_area_mm2;
label.overlap.files={};

on_target.hotspot.label_percentile=threshold.percentile;
on_target.hotspot.label_threshold_index=threshold_idx;
on_target.hotspot.label=label;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_files=etc_local_write_vertex_labels_by_hemi(decimation, vertex_index1, vertex_value, label_prefix)

if(exist('inverse_write_label','file')~=2)
    error('Cannot find inverse_write_label.m on the MATLAB path.');
end;

vertex_index1=vertex_index1(:);
vertex_value=vertex_value(:);
if(isempty(vertex_index1))
    label_files={};
    return;
end;
if(length(vertex_value)~=length(vertex_index1))
    error('Label vertex values must have the same length as label vertex indices.');
end;

label_files={};
label_idx=0;

for hemi_idx=1:length(decimation.hemi)
    hemi_name=lower(decimation.hemi(hemi_idx).name);
    full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
    full_count=decimation.hemi(hemi_idx).full_vertex_count;
    hemi_full_idx1=(full_offset+1):(full_offset+full_count);

    [is_hemi_roi, hemi_roi_pos]=ismember(vertex_index1, hemi_full_idx1(:));
    if(~any(is_hemi_roi))
        continue;
    end;

    vertex_index0=hemi_roi_pos(is_hemi_roi)-1;
    vertex_coord=decimation.full_vertex_coords(vertex_index1(is_hemi_roi),:);
    hemi_value=vertex_value(is_hemi_roi);

    fn=sprintf('%s-%s.label',label_prefix,hemi_name);
    inverse_write_label(vertex_index0, vertex_coord(:,1), vertex_coord(:,2), vertex_coord(:,3), hemi_value, fn);

    label_idx=label_idx+1;
    label_files{label_idx}=fn;
end;

% Also write one combined lh+rh label for renderers that show both
% hemispheres as one concatenated surface. The file intentionally ends in
% "-lh.label" because some label readers infer format from that suffix.
combined_vertex_index0=vertex_index1-1;
combined_vertex_coord=decimation.full_vertex_coords(vertex_index1,:);
fn=sprintf('%s_both-lh.label',label_prefix);
inverse_write_label(combined_vertex_index0, combined_vertex_coord(:,1), combined_vertex_coord(:,2), combined_vertex_coord(:,3), vertex_value, fn);

label_idx=label_idx+1;
label_files{label_idx}=fn;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hemi_count=etc_local_count_vertices_by_hemi(decimation, vertex_index1)

vertex_index1=vertex_index1(:);
hemi_count=struct([]);

for hemi_idx=1:length(decimation.hemi)
    hemi_name=lower(decimation.hemi(hemi_idx).name);
    full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
    full_count=decimation.hemi(hemi_idx).full_vertex_count;
    hemi_full_idx1=(full_offset+1):(full_offset+full_count);
    hemi_count(hemi_idx).name=hemi_name;
    hemi_count(hemi_idx).nvertices=sum(ismember(vertex_index1,hemi_full_idx1(:)));
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_data=etc_local_read_label_file(subject, label_file)

if(exist(label_file,'file')==2)
    label_data=read_label('',label_file);
else
    label_name=regexprep(label_file,'\.label$','');
    label_data=read_label(subject,label_name);
end;

if(isempty(label_data))
    error('Could not read label [%s].',label_file);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hemi=etc_local_infer_hemi_from_label_file(label_file)

[~,name,~]=fileparts(label_file);
name=lower(name);
hemi='';

if(~isempty(regexp(name,'(^|[._-])lh([._-]|$)','once')))
    hemi='lh';
elseif(~isempty(regexp(name,'(^|[._-])rh([._-]|$)','once')))
    hemi='rh';
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx=etc_local_full_indices_for_hemi(decimation, hemi, target_coord)

hemi=lower(hemi);
idx=[];

if(strcmpi(hemi,'auto'))
    hemi_dist=zeros(length(decimation.hemi),1);
    for hemi_idx=1:length(decimation.hemi)
        tmp_idx=(decimation.hemi(hemi_idx).full_vertex_offset+1):(decimation.hemi(hemi_idx).full_vertex_offset+decimation.hemi(hemi_idx).full_vertex_count);
        tmp_coord=decimation.full_vertex_coords(tmp_idx,:);
        hemi_dist(hemi_idx)=min(sum((tmp_coord-repmat(target_coord(:).',length(tmp_idx),1)).^2,2));
    end;
    [~,hemi_idx]=min(hemi_dist);
    idx=(decimation.hemi(hemi_idx).full_vertex_offset+1):(decimation.hemi(hemi_idx).full_vertex_offset+decimation.hemi(hemi_idx).full_vertex_count);
    idx=idx(:);
    return;
end;

for hemi_idx=1:length(decimation.hemi)
    if(strcmpi(hemi,'both') || strcmpi(decimation.hemi(hemi_idx).name,hemi))
        tmp_idx=(decimation.hemi(hemi_idx).full_vertex_offset+1):(decimation.hemi(hemi_idx).full_vertex_offset+decimation.hemi(hemi_idx).full_vertex_count);
        idx=cat(1,idx(:),tmp_idx(:));
    end;
end;

if(isempty(idx))
    error('Unsupported on-target hemisphere [%s]. Use ''both'', ''auto'', ''lh'', or ''rh''.',hemi);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function summary=etc_local_value_summary(value)

value=value(:);
value=value(isfinite(value));

summary=struct;
summary.n=length(value);
if(isempty(value))
    summary.mean=nan;
    summary.median=nan;
    summary.std=nan;
    summary.min=nan;
    summary.max=nan;
    summary.p90=nan;
    summary.p95=nan;
    summary.p99=nan;
    return;
end;

value_sorted=sort(value);
summary.mean=mean(value);
summary.median=median(value);
summary.std=std(value);
summary.min=min(value);
summary.max=max(value);
summary.p90=etc_local_sorted_percentile(value_sorted,90);
summary.p95=etc_local_sorted_percentile(value_sorted,95);
summary.p99=etc_local_sorted_percentile(value_sorted,99);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view_state=etc_local_render_efield(efield_to_render, output_stem, overlay_threshold, view_state)

global etc_render_fsbrain

if(nargin<3)
    overlay_threshold=[];
end;
if(nargin<4)
    view_state=[];
end;
flag_init_view=isempty(view_state);

[overlay_vertex, overlay_value]=etc_local_get_surface_overlay(efield_to_render);
overlay_value=etc_local_sanitize_overlay_value(overlay_value,'surface');
overlay_threshold=etc_local_resolve_overlay_threshold(overlay_value, overlay_threshold);

etc_render_fsbrain.overlay_vertex=overlay_vertex;
etc_render_fsbrain.overlay_value=overlay_value;

etc_render_fsbrain.overlay_smooth=5;
etc_render_fsbrain.overlay_source=1;
etc_render_fsbrain.overlay_threshold=overlay_threshold;
etc_render_fsbrain.overlay_value_flag_pos=1;
etc_render_fsbrain.overlay_value_flag_neg=0;
etc_render_fsbrain.overlay_exclude=[];
etc_render_fsbrain.overlay_include=[];
etc_render_fsbrain.overlay_fixval_flag=0;
etc_render_fsbrain.overlay_regrid_flag=0;
etc_render_fsbrain.overlay_regrid_zero_flag=0;
etc_render_fsbrain.flag_overlay_truncate_pos=0;
etc_render_fsbrain.flag_overlay_truncate_neg=0;
etc_render_fsbrain.overlay_Ds=[];
if(isfield(etc_render_fsbrain,'overlay_D'))
    etc_render_fsbrain.overlay_D=[];
end;

% Let the original initializer rebuild vol_A/overlay_vol_stc. Assigning
% overlay_vol_stc directly can leave the orthogonal-slice source indexing stale.
etc_render_fsbrain.overlay_vol_stc=[];
etc_render_fsbrain_overlay_vol_init;

etc_render_fsbrain.overlay_flag_render=1;
if(~flag_init_view)
    etc_local_restore_view_state(view_state);
end;
etc_render_fsbrain_handle('update_overlay_vol');
etc_render_fsbrain_handle('redraw');
etc_local_force_surface_overlay(overlay_value, overlay_threshold);
if(~flag_init_view)
    etc_local_restore_view_state(view_state);
end;
etc_render_fsbrain_handle('draw_pointer');
if(flag_init_view)
    etc_local_axis_tight_vis3d;
    view_state=etc_local_capture_view_state;
else
    etc_local_restore_view_state(view_state);
    etc_local_axis_vis3d;
end;

hgexport(gcf,sprintf('%s_init.png',output_stem), hgexport('factorystyle'),'Format','png');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_axis_tight_vis3d

global etc_render_fsbrain

if(~isfield(etc_render_fsbrain,'brain_axis'))
    return;
end;
if(~ishandle(etc_render_fsbrain.brain_axis))
    return;
end;

axis(etc_render_fsbrain.brain_axis,'tight');
axis(etc_render_fsbrain.brain_axis,'vis3d');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_axis_vis3d

global etc_render_fsbrain

if(~isfield(etc_render_fsbrain,'brain_axis'))
    return;
end;
if(~ishandle(etc_render_fsbrain.brain_axis))
    return;
end;

axis(etc_render_fsbrain.brain_axis,'vis3d');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view_state=etc_local_capture_view_state

global etc_render_fsbrain

view_state=[];
if(~isfield(etc_render_fsbrain,'brain_axis'))
    return;
end;
if(~ishandle(etc_render_fsbrain.brain_axis))
    return;
end;

ax=etc_render_fsbrain.brain_axis;
view_state.XLim=get(ax,'XLim');
view_state.YLim=get(ax,'YLim');
view_state.ZLim=get(ax,'ZLim');
view_state.CameraPosition=get(ax,'CameraPosition');
view_state.CameraTarget=get(ax,'CameraTarget');
view_state.CameraUpVector=get(ax,'CameraUpVector');
view_state.CameraViewAngle=get(ax,'CameraViewAngle');
view_state.DataAspectRatio=get(ax,'DataAspectRatio');
view_state.PlotBoxAspectRatio=get(ax,'PlotBoxAspectRatio');
view_state.View=get(ax,'View');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_restore_view_state(view_state)

global etc_render_fsbrain

if(isempty(view_state))
    return;
end;
if(~isfield(etc_render_fsbrain,'brain_axis'))
    return;
end;
if(~ishandle(etc_render_fsbrain.brain_axis))
    return;
end;

ax=etc_render_fsbrain.brain_axis;
try
    set(ax, ...
        'XLim',view_state.XLim, ...
        'YLim',view_state.YLim, ...
        'ZLim',view_state.ZLim, ...
        'DataAspectRatio',view_state.DataAspectRatio, ...
        'PlotBoxAspectRatio',view_state.PlotBoxAspectRatio, ...
        'View',view_state.View);
    set(ax, ...
        'CameraPosition',view_state.CameraPosition, ...
        'CameraTarget',view_state.CameraTarget, ...
        'CameraUpVector',view_state.CameraUpVector, ...
        'CameraViewAngle',view_state.CameraViewAngle);
    set(ax, ...
        'XLimMode','manual', ...
        'YLimMode','manual', ...
        'ZLimMode','manual', ...
        'DataAspectRatioMode','manual', ...
        'PlotBoxAspectRatioMode','manual', ...
        'CameraPositionMode','manual', ...
        'CameraTargetMode','manual', ...
        'CameraUpVectorMode','manual', ...
        'CameraViewAngleMode','manual');
catch ME
    fprintf('Could not restore fixed render view. %s\n',ME.message);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overlay_value=etc_local_sanitize_overlay_value(overlay_value, label)

bad_idx=find(~isfinite(overlay_value(:)));
if(~isempty(bad_idx))
    fprintf('%s overlay has %d non-finite values; setting them to zero before rendering.\n',label,length(bad_idx));
    overlay_value(bad_idx)=0;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [overlay_vertex, overlay_value]=etc_local_get_surface_overlay(efield_to_render)

global etc_render_fsbrain

n_render=size(etc_render_fsbrain.vertex_coords,1);
full_value=efield_to_render.Etotal(:);
overlay_vertex=[0:n_render-1].';

if(n_render==length(full_value))
    overlay_value=full_value;
    return;
end;

if(isfield(efield_to_render,'decimation'))
    decimation=efield_to_render.decimation;
    if(isfield(etc_render_fsbrain,'hemi'))
        if(ischar(etc_render_fsbrain.hemi))
            hemi_idx=find(strcmpi({decimation.hemi.name},etc_render_fsbrain.hemi),1);
            if(~isempty(hemi_idx))
                full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
                full_count=decimation.hemi(hemi_idx).full_vertex_count;
                if(n_render==full_count)
                    overlay_value=full_value((full_offset+1):(full_offset+full_count));
                    return;
                end;
            end;
        end;
    end;

    if(isfield(decimation,'full_vertex_coords'))
        idx=etc_local_nearest_vertex_index(decimation.full_vertex_coords, etc_render_fsbrain.vertex_coords, 2000);
        overlay_value=full_value(idx);
        return;
    end;
end;

error('Cannot match E-field values (%d) to active surface vertices (%d).',length(full_value),n_render);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overlay_threshold=etc_local_resolve_overlay_threshold(overlay_value, overlay_threshold)

if(isempty(overlay_threshold))
    overlay_threshold=etc_local_auto_overlay_threshold(overlay_value);
    return;
end;

flag_is_string=0;
if(exist('isstring','builtin')==5 || exist('isstring','file')==2)
    flag_is_string=isstring(overlay_threshold);
end;
if(ischar(overlay_threshold) || flag_is_string)
    mode=lower(char(overlay_threshold));
    switch(mode)
        case {'tans','tans_percentile','hotspot_percentile'}
            overlay_threshold=etc_local_percentile_overlay_threshold(overlay_value,[99 99.9]);
        otherwise
            error('Unsupported overlay_threshold mode [%s]. Use ''tans'', [], or a numeric [low high] threshold.',mode);
    end;
    return;
end;

if(~isnumeric(overlay_threshold) || numel(overlay_threshold)~=2)
    error('overlay_threshold must be [], ''tans'', or a numeric [low high] threshold.');
end;

overlay_threshold=double(overlay_threshold(:)).';
if(any(~isfinite(overlay_threshold)))
    error('overlay_threshold must contain finite values.');
end;
if(overlay_threshold(2)<=overlay_threshold(1))
    overlay_threshold=sort(overlay_threshold);
end;
if(overlay_threshold(2)<=overlay_threshold(1))
    overlay_threshold(2)=overlay_threshold(1)+eps;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overlay_threshold=etc_local_percentile_overlay_threshold(overlay_value, percentile_range)

val=overlay_value(:);
val=val(isfinite(val));
val=val(val>0);

if(isempty(val))
    overlay_threshold=[0 1];
    return;
end;

lo=etc_local_percentile(val,min(percentile_range));
hi=etc_local_percentile(val,max(percentile_range));

if(~isfinite(lo) || ~isfinite(hi) || hi<=lo)
    val=sort(val);
    lo=etc_local_sorted_percentile(val,min(percentile_range));
    hi=etc_local_sorted_percentile(val,max(percentile_range));
end;
if(hi<=lo)
    hi=max(val);
    lo=min(val);
end;
if(hi<=lo)
    hi=lo+eps;
end;

overlay_threshold=[lo hi];
fprintf('TANS overlay threshold from Etotal percentiles [%1.1f %1.1f] = [%1.3g %1.3g]\n', ...
    min(percentile_range), max(percentile_range), overlay_threshold(1), overlay_threshold(2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overlay_threshold=etc_local_auto_overlay_threshold(overlay_value)

val=overlay_value(:);
val=val(isfinite(val));
val=val(val>0);

if(isempty(val))
    overlay_threshold=[0 1];
    return;
end;

val=sort(val);
lo=etc_local_sorted_percentile(val,50);
hi=etc_local_sorted_percentile(val,95);

if(hi<=lo)
    hi=max(val);
    lo=min(val);
end;
if(hi<=lo)
    hi=lo+eps;
end;

overlay_threshold=[lo hi];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value=etc_local_sorted_percentile(sorted_value, pct)

idx=round(length(sorted_value).*pct./100);
idx=max(1,min(length(sorted_value),idx));
value=sorted_value(idx);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value=etc_local_percentile(value_in, pct)

value_in=value_in(:);
value_in=value_in(isfinite(value_in));
if(isempty(value_in))
    value=nan;
    return;
end;

if(exist('prctile','file')==2)
    value=prctile(value_in,pct);
else
    value=etc_local_sorted_percentile(sort(value_in),pct);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_force_surface_overlay(overlay_value, overlay_threshold)

global etc_render_fsbrain

ovs=overlay_value(:);

if(~isempty(etc_render_fsbrain.overlay_smooth))
    try
        [ovs,~,~,etc_render_fsbrain.overlay_Ds]=inverse_smooth('', ...
            'vertex',etc_render_fsbrain.vertex_coords', ...
            'face',etc_render_fsbrain.faces', ...
            'value_idx',[1:length(overlay_value)].', ...
            'value',overlay_value(:), ...
            'step',etc_render_fsbrain.overlay_smooth, ...
            'flag_fixval',0, ...
            'exc_vertex',[], ...
            'inc_vertex',[], ...
            'flag_regrid',0, ...
            'flag_regrid_zero',0, ...
            'Ds',[], ...
            'n_ratio',1);
    catch ME
        fprintf('surface smoothing failed; using unsmoothed surface overlay. %s\n',ME.message);
        ovs=overlay_value(:);
    end;
end;

if(isfield(etc_render_fsbrain,'fvdata'))
    fvdata=etc_render_fsbrain.fvdata;
else
    fvdata=repmat(etc_render_fsbrain.default_solid_color,[size(etc_render_fsbrain.vertex_coords,1),1]);
end;
if(size(fvdata,1)~=length(ovs))
    fvdata=repmat(etc_render_fsbrain.default_solid_color,[length(ovs),1]);
end;

c_idx=find(ovs(:)>=min(overlay_threshold));
if(~isempty(c_idx))
    fvdata(c_idx,:)=inverse_get_color(etc_render_fsbrain.overlay_cmap,ovs(c_idx),max(overlay_threshold),min(overlay_threshold));
end;

etc_render_fsbrain.fvdata=fvdata;
etc_render_fsbrain.ovs=ovs;

if(isfield(etc_render_fsbrain,'h'))
    if(ishandle(etc_render_fsbrain.h))
        set(etc_render_fsbrain.h,'FaceVertexCData',fvdata,'CDataMapping','direct','facecolor','interp');
        drawnow;
    end;
end;

fprintf('surface overlay: n=%d, threshold=[%1.3g %1.3g], max=%1.3g\n',length(ovs),min(overlay_threshold),max(overlay_threshold),max(ovs));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [status, bem_obj, brain_decimation]=etc_local_prepare_or_load_decimated_bem(subjects_dir, subject, file_surf, output_file_surf, tissue_name, tissue_conductivity, tissue_enclosing, file_bem, path_bem, surf_category, surface_reduce_fraction, flag_force_rebuild, nearest_chunk_size)

status=0;
info_file=fullfile(path_bem,'decimation_info.mat');

if(~exist(path_bem,'dir'))
    mkdir(path_bem);
end;

if(~flag_force_rebuild && etc_local_decimated_bem_exists(path_bem, subject, output_file_surf, file_bem, info_file))
    load(info_file,'brain_decimation');
    brain_decimation.flag_decimated_surfaces_rebuilt=0;
    bem_obj=etc_local_load_bem_obj(path_bem, subject, output_file_surf, surf_category, brain_decimation);
    status=1;
    fprintf('Using existing decimated BEM surfaces in [%s].\n',path_bem);
    return;
end;

bem_obj=struct([]);
surface_decimation=struct([]);
brain_decimation=struct;
brain_decimation.reduce_fraction_by_surface=surface_reduce_fraction(:);
brain_decimation.hemi=struct([]);
brain_decimation.full_vertex_coords=[];
brain_decimation.full_faces=[];
brain_decimation.sparse_vertex_coords=[];
brain_decimation.sparse_faces=[];
brain_decimation.sparse_full_vertex_index=[];
brain_decimation.total_full_bem_faces=0;
brain_decimation.total_sparse_bem_faces=0;
brain_decimation.flag_decimated_surfaces_rebuilt=1;

for f_idx=1:length(file_surf)
    [vertex, face0, source_files]=etc_local_read_source_surface(subjects_dir, subject, file_surf{f_idx});
    [vertex_sparse, face_sparse1]=etc_local_reduce_surface(vertex, face0, surface_reduce_fraction(f_idx));

    TR=triangulation(face_sparse1,vertex_sparse);
    surf_center=incenter(TR);
    surf_norm=faceNormal(TR);

    P=vertex_sparse;
    t=face_sparse1;
    normals=surf_norm;
    save(fullfile(path_bem,sprintf('%s_%s.mat',subject,output_file_surf{f_idx})),'P','t','normals');

    bem_obj(f_idx).filename=source_files{end};
    bem_obj(f_idx).vertex=vertex_sparse;
    bem_obj(f_idx).face=face_sparse1;
    bem_obj(f_idx).surf_center=surf_center;
    bem_obj(f_idx).surf_norm=surf_norm;
    bem_obj(f_idx).filemat=sprintf('%s_%s.mat',subject,output_file_surf{f_idx});
    if(~isempty(surf_category))
        bem_obj(f_idx).category=surf_category{f_idx};
    end;

    surface_decimation(f_idx).name=output_file_surf{f_idx};
    surface_decimation(f_idx).source_files=source_files;
    surface_decimation(f_idx).reduce_fraction=surface_reduce_fraction(f_idx);
    surface_decimation(f_idx).full_vertex_count=size(vertex,1);
    surface_decimation(f_idx).full_face_count=size(face0,1);
    surface_decimation(f_idx).sparse_vertex_count=size(vertex_sparse,1);
    surface_decimation(f_idx).sparse_face_count=size(face_sparse1,1);

    brain_decimation.total_full_bem_faces=brain_decimation.total_full_bem_faces+size(face0,1);
    brain_decimation.total_sparse_bem_faces=brain_decimation.total_sparse_bem_faces+size(face_sparse1,1);

    if(strcmpi(output_file_surf{f_idx},'gm_lh') || strcmpi(output_file_surf{f_idx},'gm_rh'))
        full_vertex_offset=size(brain_decimation.full_vertex_coords,1);
        sparse_vertex_offset=size(brain_decimation.sparse_vertex_coords,1);

        sparse_to_full_local=etc_local_sparse_to_full_index(vertex, vertex_sparse, nearest_chunk_size);
        sparse_to_full_global=sparse_to_full_local(:)+full_vertex_offset;

        brain_decimation.full_faces=cat(1,brain_decimation.full_faces,face0+full_vertex_offset);
        brain_decimation.full_vertex_coords=cat(1,brain_decimation.full_vertex_coords,vertex);
        brain_decimation.sparse_faces=cat(1,brain_decimation.sparse_faces,face_sparse1-1+sparse_vertex_offset);
        brain_decimation.sparse_vertex_coords=cat(1,brain_decimation.sparse_vertex_coords,vertex_sparse);
        brain_decimation.sparse_full_vertex_index=cat(1,brain_decimation.sparse_full_vertex_index,sparse_to_full_global);

        hemi_idx=length(brain_decimation.hemi)+1;
        brain_decimation.hemi(hemi_idx).name=regexprep(output_file_surf{f_idx},'^gm_','');
        brain_decimation.hemi(hemi_idx).source_files=source_files;
        brain_decimation.hemi(hemi_idx).full_vertex_offset=full_vertex_offset;
        brain_decimation.hemi(hemi_idx).full_vertex_count=size(vertex,1);
        brain_decimation.hemi(hemi_idx).full_face_count=size(face0,1);
        brain_decimation.hemi(hemi_idx).sparse_vertex_offset=sparse_vertex_offset;
        brain_decimation.hemi(hemi_idx).sparse_vertex_count=size(vertex_sparse,1);
        brain_decimation.hemi(hemi_idx).sparse_face_count=size(face_sparse1,1);
        brain_decimation.hemi(hemi_idx).sparse_full_vertex_index=sparse_to_full_global;
    end;
end;

brain_decimation.full_nvertices=size(brain_decimation.full_vertex_coords,1);
brain_decimation.full_nfaces=size(brain_decimation.full_faces,1);
brain_decimation.sparse_nvertices=size(brain_decimation.sparse_vertex_coords,1);
brain_decimation.sparse_nfaces=size(brain_decimation.sparse_faces,1);
brain_decimation.surface=surface_decimation;

etc_local_write_tissue_index(path_bem, file_bem, subject, output_file_surf, tissue_name, tissue_conductivity, tissue_enclosing);
save(info_file,'brain_decimation','surface_decimation');

status=1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag=etc_local_decimated_bem_exists(path_bem, subject, output_file_surf, file_bem, info_file)

flag=exist(info_file,'file')==2 && exist(fullfile(path_bem,file_bem),'file')==2;
for idx=1:length(output_file_surf)
    flag=flag && exist(fullfile(path_bem,sprintf('%s_%s.mat',subject,output_file_surf{idx})),'file')==2;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bem_obj=etc_local_load_bem_obj(path_bem, subject, output_file_surf, surf_category, brain_decimation)

bem_obj=struct([]);
for idx=1:length(output_file_surf)
    S=load(fullfile(path_bem,sprintf('%s_%s.mat',subject,output_file_surf{idx})));
    TR=triangulation(S.t,S.P);

    bem_obj(idx).filename=brain_decimation.surface(idx).source_files{end};
    bem_obj(idx).vertex=S.P;
    bem_obj(idx).face=S.t;
    bem_obj(idx).surf_center=incenter(TR);
    if(isfield(S,'normals'))
        bem_obj(idx).surf_norm=S.normals;
    else
        bem_obj(idx).surf_norm=faceNormal(TR);
    end;
    bem_obj(idx).filemat=sprintf('%s_%s.mat',subject,output_file_surf{idx});
    if(~isempty(surf_category))
        bem_obj(idx).category=surf_category{idx};
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [status, t, P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC]=etc_local_load_or_prepare_model(file_bem,path_bem,flag_force_rebuild)

status=0;
file_mesh='CombinedMesh.mat';
file_meshp='CombinedMeshP.mat';
fn_mesh=fullfile(path_bem,file_mesh);
fn_meshp=fullfile(path_bem,file_meshp);

if(~flag_force_rebuild && exist(fn_mesh,'file')==2 && exist(fn_meshp,'file')==2)
    fprintf('Loading prepared decimated BEM model from [%s].\n',path_bem);
    S=load(fn_mesh);
    SP=load(fn_meshp);
    t=S.t;
    P=S.P;
    normals=S.normals;
    Center=S.Center;
    Area=S.Area;
    Indicator=S.Indicator;
    name=S.name;
    tissue=S.tissue;
    cond=S.cond;
    enclosingTissueIdx=S.enclosingTissueIdx;
    condin=S.condin;
    condout=S.condout;
    contrast=S.contrast;
    tneighbor=SP.tneighbor;
    RnumberE=SP.RnumberE;
    ineighborE=SP.ineighborE;
    EC=SP.EC;
    status=1;
    return;
end;

pdir=pwd;
cd(path_bem);
[status, t, P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC] =etc_tms_efield_prep_model(file_bem,'path_tissue_mesh',path_bem,'file_mesh',file_mesh,'file_meshp',file_meshp);
cd(pdir);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertex, face0, source_files]=etc_local_read_source_surface(subjects_dir, subject, file_entry)

vertex=[];
face0=[];
source_files=cell(size(file_entry));

for entry_idx=1:length(file_entry)
    if(length(file_entry)==1)
        source_file=fullfile(subjects_dir,subject,'bem',file_entry{entry_idx});
    else
        source_file=fullfile(subjects_dir,subject,file_entry{entry_idx});
    end;
    source_file=etc_local_resolve_subject_surface_file(subjects_dir, subject, source_file, file_entry{entry_idx});

    [v_tmp,f_tmp]=read_surf(source_file);
    source_files{entry_idx}=source_file;
    face0=cat(1,face0,f_tmp+size(vertex,1));
    vertex=cat(1,vertex,v_tmp);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source_file=etc_local_resolve_subject_surface_file(subjects_dir, subject, source_file, file_name)

if(exist(source_file,'file')==2)
    return;
end;

switch(file_name)
    case 'outer_skin.surf'
        fallback=fullfile(subjects_dir,subject,'bem','watershed',sprintf('%s_outer_skin_surface',subject));
    case 'outer_skull.surf'
        fallback=fullfile(subjects_dir,subject,'bem','watershed',sprintf('%s_outer_skull_surface',subject));
    case 'inner_skull.surf'
        fallback=fullfile(subjects_dir,subject,'bem','watershed',sprintf('%s_inner_skull_surface',subject));
    otherwise
        fallback='';
end;

if(~isempty(fallback) && exist(fallback,'file')==2)
    source_file=fallback;
    return;
end;

error('Cannot find source surface [%s].',source_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertex_sparse, face_sparse1]=etc_local_reduce_surface(vertex, face0, reduce_fraction)

if(reduce_fraction>=1)
    vertex_sparse=vertex;
    face_sparse1=face0+1;
else
    [face_sparse1, vertex_sparse]=reducepatch(double(face0)+1, double(vertex), reduce_fraction);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etc_local_write_tissue_index(path_bem, file_bem, subject, output_file_surf, tissue_name, tissue_conductivity, tissue_enclosing)

fp=fopen(fullfile(path_bem,file_bem),'w');
if(fp<0)
    error('Cannot write tissue index [%s].',fullfile(path_bem,file_bem));
end;

for tissue_idx=1:length(tissue_name)
    fprintf(fp,'>%s : %s_%s.mat : %1.4f : %s \n',tissue_name{tissue_idx},subject,output_file_surf{tissue_idx},tissue_conductivity(tissue_idx),tissue_enclosing{tissue_idx});
end;

fclose(fp);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx=etc_local_sparse_to_full_index(full_vertices, sparse_vertices, nearest_chunk_size)

[is_exact, idx_exact]=ismember(sparse_vertices, full_vertices, 'rows');
idx=zeros(size(sparse_vertices,1),1);
idx(is_exact)=idx_exact(is_exact);

if(any(~is_exact))
    idx_missing=etc_local_nearest_vertex_index(full_vertices, sparse_vertices(~is_exact,:), nearest_chunk_size);
    idx(~is_exact)=idx_missing;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx=etc_local_nearest_vertex_index(reference_vertices, query_vertices, chunk_size)

reference_vertices=double(reference_vertices);
query_vertices=double(query_vertices);

if(exist('knnsearch','file')==2)
    idx=knnsearch(reference_vertices, query_vertices);
    idx=idx(:);
    return;
end;

idx=zeros(size(query_vertices,1),1);
reference_norm=sum(reference_vertices.^2,2).';

for q_idx=1:chunk_size:size(query_vertices,1)
    q_end=min(q_idx+chunk_size-1,size(query_vertices,1));
    q=query_vertices(q_idx:q_end,:);
    d2=bsxfun(@plus,sum(q.^2,2),reference_norm)-2.*(q*reference_vertices.');
    [~,idx(q_idx:q_end)]=min(d2,[],2);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tissue_to_plot=etc_local_hemi_to_tissue(hemi)

if(iscell(hemi))
    tissue_to_plot=cell(size(hemi));
    for hemi_idx=1:length(hemi)
        tissue_to_plot{hemi_idx}=etc_local_one_hemi_to_tissue(hemi{hemi_idx});
    end;
else
    tissue_to_plot=etc_local_one_hemi_to_tissue(hemi);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tissue_name=etc_local_one_hemi_to_tissue(hemi)

switch(lower(hemi))
    case 'lh'
        tissue_name='GM_LH';
    case 'rh'
        tissue_name='GM_RH';
    otherwise
        error('Unsupported hemisphere [%s]. Use ''lh'' or ''rh''.', hemi);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decimation=etc_local_add_full_to_sparse_map(decimation, nearest_chunk_size, smooth_step)

if(nargin<3 || isempty(smooth_step))
    smooth_step=5;
end;

full_to_sparse_row_index=zeros(decimation.full_nvertices,1);

for hemi_idx=1:length(decimation.hemi)
    full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
    full_count=decimation.hemi(hemi_idx).full_vertex_count;
    sparse_offset=decimation.hemi(hemi_idx).sparse_vertex_offset;
    sparse_count=decimation.hemi(hemi_idx).sparse_vertex_count;

    full_idx=(full_offset+1):(full_offset+full_count);
    sparse_idx=(sparse_offset+1):(sparse_offset+sparse_count);

    nearest_local=etc_local_nearest_vertex_index(decimation.sparse_vertex_coords(sparse_idx,:), decimation.full_vertex_coords(full_idx,:), nearest_chunk_size);
    full_to_sparse_row_index(full_idx)=nearest_local(:)+sparse_offset;
end;

decimation.full_to_sparse_row_index=full_to_sparse_row_index;
decimation.full_to_sparse_row_index_method='nearest';
decimation.full_interpolation_method='nearest_inverse_smooth';
decimation.full_interpolation_smooth_step=smooth_step;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stc_info=etc_local_write_efield_stc(efield_full, decimation, output_dir, output_stem)

if(exist('inverse_write_stc','file')~=2)
    error('inverse_write_stc.m is required to write STC files but is not on the MATLAB path.');
end;
if(~isfield(decimation,'hemi') || isempty(decimation.hemi))
    error('brain_decimation.hemi is required to split full E-field values into lh/rh STC files.');
end;
if(~isfield(decimation,'full_nvertices'))
    error('brain_decimation.full_nvertices is required for STC export.');
end;

scalar_names={'Etotal','Enormal','Etangent'};
stc_info=struct;
stc_info.output_dir=output_dir;
stc_info.output_stem=output_stem;
stc_info.scalar_names=scalar_names;
stc_info.epoch_begin_latency=0;
stc_info.sample_period=1;
stc_info.files=struct;

fprintf('Writing native STC E-field maps [%s]...\n',output_stem);
for scalar_idx=1:length(scalar_names)
    scalar_name=scalar_names{scalar_idx};
    value=etc_local_get_efield_stc_value(efield_full, scalar_name);
    if(length(value)~=decimation.full_nvertices)
        error('%s length %d does not match full_nvertices %d.', scalar_name, length(value), decimation.full_nvertices);
    end;

    finite_value=value(isfinite(value));
    stc_info.files.(scalar_name)=struct;
    stc_info.files.(scalar_name).min=min(finite_value);
    stc_info.files.(scalar_name).max=max(finite_value);

    for hemi_idx=1:length(decimation.hemi)
        hemi_name=lower(decimation.hemi(hemi_idx).name);
        if(~strcmp(hemi_name,'lh') && ~strcmp(hemi_name,'rh'))
            continue;
        end;
        full_offset=decimation.hemi(hemi_idx).full_vertex_offset;
        full_count=decimation.hemi(hemi_idx).full_vertex_count;
        full_idx=[full_offset+1:full_offset+full_count].';
        local_vertex0=int32([0:full_count-1].');
        stc_value=double(value(full_idx));
        stc_file=fullfile(output_dir,sprintf('%s_%s-%s.stc',output_stem,scalar_name,hemi_name));

        inverse_write_stc(stc_value(:), local_vertex0, 0, 1, stc_file);

        stc_info.files.(scalar_name).(hemi_name)=stc_file;
        stc_info.files.(scalar_name).(sprintf('n_%s_vertices',hemi_name))=full_count;
        fprintf('  %s-%s: n=%d, range=[%1.3f %1.3f]\n',scalar_name,hemi_name,full_count,min(stc_value),max(stc_value));
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value=etc_local_get_efield_stc_value(efield_full, scalar_name)

switch(lower(scalar_name))
    case 'etotal'
        if(isfield(efield_full,'Etotal'))
            value=efield_full.Etotal;
        elseif(isfield(efield_full,'E'))
            value=sqrt(sum(double(efield_full.E).^2,2));
        else
            error('efield does not contain Etotal or vector E.');
        end;
    case 'enormal'
        if(isfield(efield_full,'Enormal'))
            value=efield_full.Enormal;
        else
            error('efield does not contain Enormal.');
        end;
    case 'etangent'
        if(isfield(efield_full,'Etangent'))
            value=efield_full.Etangent;
        elseif(isfield(efield_full,'Etotal') && isfield(efield_full,'Enormal'))
            value=sqrt(max(double(efield_full.Etotal).^2-double(efield_full.Enormal).^2,0));
        else
            error('efield does not contain Etangent, or Etotal plus Enormal.');
        end;
    otherwise
        error('Unsupported STC scalar [%s].',scalar_name);
end;

if(size(value,2)>1)
    value=sqrt(sum(double(value).^2,2));
else
    value=double(value(:));
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [efield_status,efield_sparse]=etc_local_run_efield_solve(bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords, tissue_to_plot, strcoil)

[efield_status,efield_sparse]=etc_tms_efield_surf( ...
    bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, ...
    enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ...
    ineighborE, EC, coords, 'tissue_to_plot', tissue_to_plot, 'strcoil', strcoil);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xfm_m,xfm_mm,coil_center,coil_orientation]=etc_local_scalp_coil_xfm_goto(scalp_center,scalp_normal,scalp_up,tms_coil_xfm_reference)

coil_center=double(scalp_center(:)).';
outward_normal=double(scalp_normal(:)).';
outward_normal=outward_normal./norm(outward_normal);
coil_orientation=-outward_normal;
coil_axis=coil_orientation(:);

coil_up=double(scalp_up(:));
coil_up=coil_up-sum(coil_up.*coil_axis).*coil_axis;
if(norm(coil_up)<eps)
    error('Scalp candidate has no valid tangent handle direction.');
end;
coil_up=coil_up./norm(coil_up);
coil_x=cross(coil_up, -coil_axis);
coil_x=coil_x./norm(coil_x);

xfm_m=eye(4);
xfm_m(1:3,1)=coil_x;
xfm_m(1:3,2)=coil_up;
xfm_m(1:3,3)=-coil_axis;
xfm_m(1:3,4)=coil_center(:)./1e3;
xfm_mm=etc_local_tms_xfm_m_to_mm(xfm_m);

% Keep the reference argument in the interface so this helper can replace
% target-derived placement without changing downstream transform handling.
if(nargin>=4 && isempty(tms_coil_xfm_reference))
    error('A reference TMS coil transform is required.');
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function efield=etc_local_synthesize_efield_from_basis(efield0, efield90, rotation_deg)

if(~isfield(efield0,'E') || ~isfield(efield90,'E'))
    error('Basis E-field structures must contain vector field E.');
end;
if(size(efield0.E,1)~=size(efield90.E,1) || size(efield0.E,2)~=size(efield90.E,2))
    error('Basis E-fields have incompatible dimensions.');
end;

c=cosd(rotation_deg);
s=sind(rotation_deg);
efield=efield0;
efield.E=c.*double(efield0.E)+s.*double(efield90.E);
efield.Etotal=sqrt(sum(efield.E.^2,2));

if(isfield(efield0,'Enormal') && isfield(efield90,'Enormal'))
    efield.Enormal=c.*double(efield0.Enormal)+s.*double(efield90.Enormal);
else
    error('Basis E-field structures must contain Enormal for TANS synthesis.');
end;
efield.Etangent=sqrt(max(efield.Etotal.^2-efield.Enormal.^2,0));
efield.coil_rotation_deg=rotation_deg;
efield.basis_synthesis=1;
efield.basis_rotation_deg=[0 90];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function row=etc_local_empty_basis_validation_row()

row=struct;
row.candidate_id=nan;
row.rotation_deg=nan;
row.vector_rms_relative_error=nan;
row.vector_max_absolute_error=nan;
row.Etotal_rms_relative_error=nan;
row.Enormal_rms_relative_error=nan;
row.exact_solve_elapsed_sec=nan;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rows=etc_local_validate_basis_candidate(candidate_id, target_coord, rotation_deg, tms_coil_xfm_goto, tms_coil_xfm, strcoil_base, bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ineighborE, EC, coords, tissue_to_plot, basis_entry)

rows=repmat(etc_local_empty_basis_validation_row(),numel(rotation_deg),1);
for rotation_idx=1:numel(rotation_deg)
    angle=rotation_deg(rotation_idx);
    xfm=etc_local_tms_rotate_about_axis(tms_coil_xfm_goto,angle);
    [strcoil_exact,~]=etc_tms_target_xfm_apply( ...
        strcoil_base, xfm(1:3,4).*1e3, -xfm(1:3,3), xfm(1:3,2), ...
        tms_coil_xfm, xfm);
    tic;
    [~,exact]=etc_local_run_efield_solve( ...
        bem_t, bem_P, normals, Center, Area, Indicator, name, tissue, cond, ...
        enclosingTissueIdx, condin, condout, contrast, tneighbor, RnumberE, ...
        ineighborE, EC, coords, tissue_to_plot, strcoil_exact);
    elapsed=toc;
    synthesized=etc_local_synthesize_efield_from_basis( ...
        basis_entry.efield0,basis_entry.efield90,angle);
    delta=double(exact.E)-double(synthesized.E);
    exact_norm=max(norm(double(exact.E(:))),eps);
    rows(rotation_idx).candidate_id=candidate_id;
    rows(rotation_idx).rotation_deg=angle;
    rows(rotation_idx).vector_rms_relative_error=sqrt(mean(delta(:).^2))/exact_norm;
    rows(rotation_idx).vector_max_absolute_error=max(abs(delta(:)));
    rows(rotation_idx).Etotal_rms_relative_error= ...
        sqrt(mean((double(exact.Etotal(:))-double(synthesized.Etotal(:))).^2))/max(norm(double(exact.Etotal(:))),eps);
    rows(rotation_idx).Enormal_rms_relative_error= ...
        sqrt(mean((double(exact.Enormal(:))-double(synthesized.Enormal(:))).^2))/max(norm(double(exact.Enormal(:))),eps);
    rows(rotation_idx).exact_solve_elapsed_sec=elapsed;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function efield_full=etc_local_expand_efield_to_full(efield_sparse, decimation, full_nvertices, smooth_step)

if(nargin<4 || isempty(smooth_step))
    smooth_step=5;
end;

efield_full=efield_sparse;
efield_full.vertices=int32([0:full_nvertices-1].');
efield_full.flag_sparse=0;
efield_full.decimation=decimation;
efield_full.decimation.full_interpolation_method='nearest_inverse_smooth';
efield_full.decimation.full_interpolation_smooth_step=smooth_step;

smooth_D=[];
if(isfield(efield_sparse,'E'))
    [efield_full.E,smooth_D]=etc_local_smooth_sparse_field_to_full(efield_sparse.E, decimation, full_nvertices, smooth_step, smooth_D);
end;
if(isfield(efield_sparse,'Enormal'))
    [efield_full.Enormal,smooth_D]=etc_local_smooth_sparse_field_to_full(efield_sparse.Enormal, decimation, full_nvertices, smooth_step, smooth_D);
end;
if(isfield(efield_full,'E'))
    efield_full.Etotal=sqrt(sum(efield_full.E.^2,2));
elseif(isfield(efield_sparse,'Etotal'))
    [efield_full.Etotal,smooth_D]=etc_local_smooth_sparse_field_to_full(efield_sparse.Etotal, decimation, full_nvertices, smooth_step, smooth_D);
end;
if(isfield(efield_full,'Etotal') && isfield(efield_full,'Enormal'))
    efield_full.Etangent=sqrt(max(efield_full.Etotal.^2-efield_full.Enormal.^2,0));
elseif(isfield(efield_sparse,'Etangent'))
    [efield_full.Etangent,smooth_D]=etc_local_smooth_sparse_field_to_full(efield_sparse.Etangent, decimation, full_nvertices, smooth_step, smooth_D);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [full_value,smooth_D]=etc_local_smooth_sparse_field_to_full(field_value, decimation, full_nvertices, smooth_step, smooth_D)

if(nargin<5)
    smooth_D=[];
end;

field_value=double(field_value);
if(size(field_value,1)==full_nvertices)
    full_value=field_value;
    return;
end;

nearest_value=field_value(decimation.full_to_sparse_row_index,:);

if(exist('inverse_smooth','file')~=2 || smooth_step<=0)
    full_value=nearest_value;
    return;
end;

try
    [full_value,~,~,~,smooth_D]=inverse_smooth('', ...
        'vertex',decimation.full_vertex_coords', ...
        'face',double(decimation.full_faces)', ...
        'value_idx',[1:full_nvertices].', ...
        'value',nearest_value, ...
        'step',smooth_step, ...
        'flag_fixval',0, ...
        'exc_vertex',[], ...
        'inc_vertex',[], ...
        'flag_regrid',0, ...
        'flag_regrid_zero',0, ...
        'flag_scale',0, ...
        'D',smooth_D, ...
        'n_ratio',1);
catch ME
    fprintf('inverse_smooth full-density interpolation failed; using nearest fill. %s\n',ME.message);
    full_value=nearest_value;
end;

end
