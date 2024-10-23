clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri_fus/subjects/'); %for MAC/Linux

%source space
file_source_fif='/Users/fhlin/workspace/eegmri_fus/subjects/EP142/bem/EP142-5-src.fif';

%surf_outer_skin='/Users/fhlin/workspace/eegmri_fus/subjects/EP142/bem/watershed/EP142_outer_skin_surface';
surf_outer_skull='/Users/fhlin/workspace/eegmri_fus/subjects/EP142/bem/watershed/EP142_outer_skull_surface';
surf_inner_skull='/Users/fhlin/workspace/eegmri_fus/subjects/EP142/bem/watershed/EP142_inner_skull_surface';
surf_outer_skin='/Users/fhlin/workspace/eegmri_fus/subjects/EP142/surf/lh.seghead.tri';

%electrodes
file_elec_loc='../EP142_231205/digitizer/RJM.DAT'   ;
n_eeg=32; %this one is hard-coded here. But it should be determined automatically from files in the future. 


eeg_output_name='eeg_fwd_prep_102224.mat';


h_mri_brain_lh=[];
h_mri_brain_rh=[];
h_mri_head=[];
h_mri_source_location=[];
h_mri_source_norm=[];
h_meg_array=[];


% flag_show_outer_skin=0;
% flag_show_outer_skull=1;
% flag_show_inner_skull=1;
flag_show_outer_skin=0;
flag_show_outer_skull=1;
flag_show_inner_skull=0;

flag_show_mri_brain=1;
flag_show_mri_head=1;
flag_show_meg_array=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[src] = mne_read_source_spaces(file_source_fif);

%[verts_osc, faces_osc] = mne_read_surface(surf_outer_skin);
 [verts_osk, faces_osk] = mne_read_surface(surf_outer_skull);
 [verts_isk, faces_isk] = mne_read_surface(surf_inner_skull);
 [verts_osc, faces_osc] = inverse_read_tri(surf_outer_skin);
verts_osc=verts_osc.*1e-3;

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
    load(eeg_output_name,'points_label','points');
    if(exist('points'))
        elec.elecpos=points.*1e3;
    else
        [ch_names, chid, x,y,z]=etc_read_datfile(file_elec_loc);
        points(:,1)=x;
        points(:,2)=y;
        points(:,3)=z;
        
        elec.elecpos=points.*10; %load electrodes; variable 'elec'
    end;
    
    if(exist('points_label'))
        elec.label=points_label;
    end;
else
    [ch_names, chid, x,y,z]=etc_read_datfile(file_elec_loc); %x, y, z in cm
    points(:,1)=x;
    points(:,2)=y;
    points(:,3)=z;
    
    points=points.*10.*1e-3; %x, y, z in m
    
    elec.elecpos=points.*1e3; %load electrodes; variable 'elec'; size in mm
    elec.label=ch_names;
    points_label=elec.label;
    %n_eeg=length(points_label);
    
    %elec.elecpos=elec.elecpos(1:n_eeg,:);
    %elec.label=elec.label(1:n_eeg);
    %points_label=points_label(1:n_eeg);
    
    %points=[];
end;

etc_render_topo('alpha',0.8,'vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_aux_point_coords',elec.elecpos./1000,'topo_aux_point_name',elec.label);
view(-60,20);

%save for EEG solution
if(~exist(eeg_output_name))
    save(eeg_output_name,'faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points','points_label','n_eeg');
else
    fprintf('found [%s]... all variables except sensor locations are replaced!\n',eeg_output_name);
    save(eeg_output_name,'-append','faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points','points_label','n_eeg');
end;

return;
