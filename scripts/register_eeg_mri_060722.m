clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri/');

%source space
file_source_fif='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/180322_PYW-5-src.fif';

surf_outer_skin='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/watershed/180322_PYW_outer_skin_surface';
surf_outer_skull='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/watershed/180322_PYW_outer_skull_surface';
surf_inner_skull='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/watershed/180322_PYW_inner_skull_surface';

%electrodes
file_elec_loc='../digitizer/PYW.DAT';
file_elec_name='../digitizer/PYW.ela';
n_eeg=31; %this one is hard-coded here. But it should be determined automatically from files in the future. 


eeg_output_name='eeg_fwd_prep_060722.mat';


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
    load(eeg_output_name,'points_label','points');
    if(exist('points'))
        elec.elecpos=points.*1e3;
    else
        [d0,d1,dx,dy,dz]=etc_read_datfile(file_elec_loc);
        elec.elecpos(:,1)=dx.*1e1;
        elec.elecpos(:,2)=dy.*1e1;
        elec.elecpos(:,3)=dz.*1e1;
        elec.label=d0;
    end;
    elec.label=points_label;
else
    [d0,d1,dx,dy,dz]=etc_read_datfile(file_elec_loc);
    elec.elecpos(:,1)=dx.*1e1;
    elec.elecpos(:,2)=dy.*1e1;
    elec.elecpos(:,3)=dz.*1e1;
    elec.label=d0;
    points_label=elec.label;

    
    points=[];
end;

etc_render_topo('alpha',0.2,'vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_aux_point_coords',elec.elecpos./1e3,'topo_aux_point_name',elec.label);
view(-60,20);

%save for EEG solution
if(~exist(eeg_output_name))
    save(eeg_output_name,'faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points_label','n_eeg');
else
    fprintf('found [%s]... all variables except sensor locations are replaced!\n',eeg_output_name);
    save(eeg_output_name,'-append','faces_isk','faces_osk','faces_osc','verts_isk','verts_osk','verts_osc','src','points_label','n_eeg');
end;

return;