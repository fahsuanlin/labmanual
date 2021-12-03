close all; clear all;

clear global etc_render_fsbrain;

subject='fsaverage';
surf='orig';
hemi='lh';

file_label={
    'lh.Temporo-Parieto-Occipital Junctio.label';
    'lh.Early Auditory Corte.label';
    'lh.Auditory Association Corte.label';
    };

roi_stem={
    'lh.Temporo-Parieto-OccipitalJunctio';
    'lh.EarlyAuditoryCorte';
    'lh.AuditoryAssociationCorte';
};

orig_subject={
%    's025';
    's026';
    's027';
%    's031';
    's032';
    's033';
%    's034';
%    's035';
%    's036';
%    's038';
    's041';
    's045';
    's046';
%    's047';
    's048';
    's050';
%    's051';
%    's054';
    's055';
    's057';
    };

root_path='/Users/fhlin/workspace/seeg';

%color of contacts
aux2_point_individual_color=[0.6 0.5 0.5];
selected_electrode_color=[0 0.45 0.75];
selected_contact_color=[0 1 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects/'); %for MAC/Linux
%setenv('SUBJECTS_DIR','D:\fhlin\Users\fhlin\workspace\seeg\subjects'); %for PC

subjects_dir=getenv('SUBJECTS_DIR');

mri=MRIread(sprintf('%s/%s/mri/orig.mgz',subjects_dir,subject));
%mri=etc_MRIread('D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\orig.mgz'); %for PC
%
xfm=etc_read_xfm('file_xfm',sprintf('%s/%s_post/tmp/register.dat',subjects_dir,subject));
%xfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036_post\tmp\register.dat'); %for PC

%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm',sprintf('%s/%s/mri/transforms/talairach.xfm',subjects_dir,subject)); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\transforms\talairach.xfm'); %for PC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load electrode
e_counter=1;
for subj_idx=1:length(orig_subject)
    for ll_idx=1:length(roi_stem)
        file_electrode=sprintf('electrode_tal_mri_%s_%s_120121.mat',roi_stem{ll_idx},orig_subject{subj_idx});
        
        load(sprintf('%s/%s/analysis/%s',root_path,orig_subject{subj_idx},file_electrode));
        E.coord=E.coord(E.contact_idx,:);
        E.contact_idx=1;
        E.n_contact=1;
        electrode(e_counter)=E;
        
        switch ll_idx
            case 1
                e_color=[0   0.4470   0.741];
            case 2
                e_color=[0.8500 0.3250 0.0980];
            case 3
                e_color=[0.9290 0.6940 0.1250];
        end;
        
        aux2_point_individual_color(e_counter,:)=e_color;
        aux2_point_coords(e_counter,:)=E.coord(1,:);
        aux2_point_name{e_counter}=sprintf('%s_%d',E(1).name,1);      
        
        e_counter=e_counter+1;
    end;
end;

%load label
if(~isempty(file_label))
    for ll_idx=1:length(file_label)
        [ll{ll_idx}]=inverse_read_label(file_label{ll_idx});
        [dumm,label_file_stem]=fileparts(file_label{ll_idx});
        
        label_vertex(ll{ll_idx}+1)=ll_idx;
        label_value(ll{ll_idx}+1)=ll_idx;
        label_ctab.struct_names{ll_idx}=label_file_stem;
        switch mod(ll_idx,5)
            case 1
                label_ctab.table(ll_idx,:)=[0*256   0.4470*256   0.741*256         0       ll_idx];
            case 2
                label_ctab.table(ll_idx,:)=[0.8500*256 0.3250*256 0.0980*256          0        ll_idx];
            case 3
                label_ctab.table(ll_idx,:)=[0.9290*256 0.6940*256 0.1250*256          0        ll_idx];
            case 4
                label_ctab.table(ll_idx,:)=[0.4940*256 0.1840*256 0.5560*256          0        ll_idx];
            case 5
                label_ctab.table(ll_idx,:)=[0.4660*256 0.6740*256 0.1880*256          0        ll_idx];
        end;
    end
    label_ctab.numEntries=length(file_label);
        
else
    ll=[];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
if(~isempty(electrode))
    [orig_vertex_coords, orig_faces] = read_surf(sprintf('%s/%s/surf/%s.orig',subjects_dir,subject,hemi));
    [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,hemi,surf));
    
    count=1;
    for e_idx=1:length(electrode)
        for c_idx=1:electrode(e_idx).n_contact
           
%             aux2_point_coords(count,:)=electrode(e_idx).coord(c_idx,:);
             
%             aux2_point_individual_color(count,:)=selected_electrode_color;
%             if(c_idx==electrode(e_idx).contact_idx)
%                 aux2_point_individual_color(count,:)=selected_contact_color;
%             end;
%             
%             aux2_point_name{count}=sprintf('%s_%d',electrode(e_idx).name, c_idx);;
%             count=count+1;
        end;
    end;
end;




%etc_render_fsbrain('overlay_cmap',[1.0000    0.4118    0.1608].*0.8,'overlay_smooth',[],'overlay_vertex',ll,'overlay_value',ones(length(ll),1),'overlay_threshold',[0.2 0.3],'curv_pos_color',[1 1 1].*0.5,'curv_neg_color',[1 1 1].*0.5,'alpha',0.3,'surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'topo_aux2_point_coords',aux2_point_coords,'electrode',electrode);
etc_render_fsbrain('curv_pos_color',[1 1 1].*0.5,'curv_neg_color',[1 1 1].*0.5,'alpha',0.3,'surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'topo_aux2_point_coords',aux2_point_coords,'electrode',electrode);
view(-90,30);

global etc_render_fsbrain;
etc_render_fsbrain.label_value=label_value;
etc_render_fsbrain.label_vertex=label_vertex;
etc_render_fsbrain.label_ctab=label_ctab;
etc_render_fsbrain.label_register=ones(length(file_label),1);

etc_render_fsbrain_handle('update_label');

hgexport(gcf,sprintf('all_label_tal_brain_120121_%s',hemi), hgexport('factorystyle'),'Format','png');

global etc_render_fsbrain;
etc_render_fsbrain.aux2_point_individual_color=aux2_point_individual_color;
% etc_render_fsbrain.selected_electrode_color=selected_electrode_color;
% etc_render_fsbrain.selected_contact_color=selected_contact_color;
etc_render_fsbrain.selected_electrode_flag=0;
etc_render_fsbrain.selected_contact_flag=0;

etc_render_fsbrain_handle('redraw');
etc_render_fsbrain_handle('update_label');


hgexport(gcf,sprintf('all_electrode_label_tal_brain_120121_%s',hemi), hgexport('factorystyle'),'Format','png');


return;
