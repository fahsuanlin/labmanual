clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects/'); %for MAC/Linux

%source space
file_source_fif='/Users/fhlin/workspace/seeg/subjects/s065/bem/s065-5-src.fif';

surf_outer_skin='/Users/fhlin/workspace/seeg/subjects/s065/bem/watershed/s065_outer_skin_surface';
surf_outer_skull='/Users/fhlin/workspace/seeg/subjects/s065/bem/watershed/s065_outer_skull_surface';
surf_inner_skull='/Users/fhlin/workspace/seeg/subjects/s065/bem/watershed/s065_inner_skull_surface';

%electrodes
file_elp='../digitizer/s065_HCH_20201215_01.elp';
file_pos='../digitizer/s065_HCH_20201215_01.pos';
n_eeg=32; %this one is hard-coded here. But it should be determined automatically from files in the future. 


eeg_output_name='eeg_fwd_prep_041721.mat';


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


if(exist(eeg_output_name))
    fprintf('loading [%s]...\n',eeg_output_name);
    load(eeg_output_name,'elec','points_label','points');
    if(exist('points'))
        elec.elecpos=points.*1e3;
    else
        %elec.elecpos=load(file_elec_loc).*10; %load electrodes; variable 'elec'
        %elec.elecpos(:,2)=elec.elecpos(:,2)-32; %rough manual alignment
        %elec.elecpos(:,3)=elec.elecpos(:,3)-50; %rough manual alignment
 
        [elp_names,elp_x,elp_y,elp_z]=etc_read_posfile(file_pos,'n_headerline',1);
         elp_points=cat(2,elp_x(:),elp_y(:),elp_z(:));
        elp_points=elp_points;

        elec.elecpos=elp_points;
        elec.label=elp_names;
        for ii=1:size(elp_points,1)
            elec.label{end+1}='';
        end;
        
        %elec.elecpos=elec.elecpos(1:n_eeg,:);
    end;
    elec.label=points_label;
else
    %elec.elecpos=load(file_elec_loc).*10; %load electrodes; variable 'elec'
    %elec.elecpos(:,2)=elec.elecpos(:,2)-32; %rough manual alignment
    %elec.elecpos(:,3)=elec.elecpos(:,3)-50; %rough manual alignment
    
    [elp_names,elp_x,elp_y,elp_z]=etc_read_posfile(file_pos,'n_headerline',1);
    elp_points=cat(2,elp_x(:),elp_y(:),elp_z(:));
    elp_points=elp_points;
    
    elec.elecpos=elp_points;
    elec.label=elp_names;
    for ii=1:size(elp_points,1)
        elec.label{end+1}='';
    end;
    
    points_label=elec.label;
    n_eeg=length(points_label);
    
    %elec.elecpos=elec.elecpos(1:n_eeg,:);
    %elec.label=elec.label(1:n_eeg);
    %points_label=points_label(1:n_eeg);
    
    points=[];
end;

etc_render_topo('alpha',0.2,'vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_aux_point_coords',elec.elecpos./1e2,'topo_aux_point_name',elec.label);
view(-60,20);

%save for EEG solution
if(~exist(eeg_output_name))
    save(eeg_output_name,'faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points_label','n_eeg');
else
    fprintf('found [%s]... all variables except sensor locations are replaced!\n',eeg_output_name);
    save(eeg_output_name,'-append','faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points_label','n_eeg');
end;

return;