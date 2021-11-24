close all; clear all;

subject='s031';
surf='orig';
hemi='rh';

file_label='aud_HCPMMP1-rh_s031.label';

file_electrode='electrode_101721_223845_s031.mat'; %electrode coordinates in post-op MRI

file_roi='electrodes_to_labels_061721.mat';
roi_idx=4; %auditory cortex; right hemisphere

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
if(~isempty(file_electrode))
    load(file_electrode);
else
    electrode=[];
end;

%load label
if(~isempty(file_label))
    [ll]=inverse_read_label(file_label);
else
    ll=[];
end;

%find electrode and contact index
if(~isempty(file_roi))
    load(file_roi);
    if(~isempty(electrode))
        for electrode_idx=1:length(electrode)
            if(strcmp(electrode(electrode_idx).name,roi(roi_idx).electrode_min_dist_electrode_name))
                break;
            end;
        end;
    end;
    contact_idx=roi(roi_idx).electrode_min_dist_electrode_contact;
else
    roi=[];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isempty(electrode))
    [orig_vertex_coords, orig_faces] = read_surf(sprintf('%s/%s/surf/%s.orig',subjects_dir,subject,hemi));
    [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,hemi,surf));
    
    count=1;
    for e_idx=1:length(electrode)
        for c_idx=1:electrode(e_idx).n_contact
            
            aux2_point_coords(count,:)=electrode(e_idx).coord(c_idx,:);
            
            if(strcmp(surf,'orig')|strcmp(surf,'smoothwm')|strcmp(surf,'pial'))
                
            else
                fprintf('surface <%s> not "orig"/"smoothwm"/"pial". Electrode contacts locations are updated to the nearest location of this surface.\n',etc_render_fsbrain.surf);
                
                tmp=aux2_point_coords(count,:);
                
                vv=orig_vertex_coords;
                dist=sqrt(sum((vv-repmat([tmp(1),tmp(2),tmp(3)],[size(vv,1),1])).^2,2));
                [min_dist,min_dist_idx]=min(dist);
                if(~isnan(min_dist))
                    aux2_point_coords(count,:)=vertex_coords(min_dist_idx,:);
                end;
            end;
            
            aux2_point_name{count}=sprintf('%s_%d',electrode(e_idx).name, c_idx);;
            count=count+1;
        end;
    end;
end;




etc_render_fsbrain('overlay_cmap',[1.0000    0.4118    0.1608].*0.8,'overlay_smooth',[],'overlay_vertex',ll,'overlay_value',ones(length(ll),1),'overlay_threshold',[0.2 0.3],'curv_pos_color',[1 1 1].*0.5,'curv_neg_color',[1 1 1].*0.5,'alpha',0.3,'surf',surf,'hemi','rh','subject','s031','vol',mri,'talxfm',(talxfm),'topo_aux2_point_coords',aux2_point_coords,'electrode',electrode,'electrode_idx',electrode_idx,'electrode_contact_idx',contact_idx);
view(90,30);

global etc_render_fsbrain;
etc_render_fsbrain.aux2_point_color=[0.6 0.5 0.5];
etc_render_fsbrain.selected_electrode_color=[0 0.45 0.75];
etc_render_fsbrain.selected_contact_color=[0 1 1];

etc_render_fsbrain_handle('redraw');
return;
