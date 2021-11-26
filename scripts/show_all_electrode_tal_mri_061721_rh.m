close all; clear all;

clear global etc_render_fsbrain;

subject='fsaverage';
surf='orig';
hemi='rh';

file_label='aud_HCPMMP1-rh.label';

roi_stem='aud_HCPMMP1-rh';

orig_subject={
    's025';
    %    's026';
    %    's027';
    's031';
    %     %    's032';
    %     %    's033';
    's034';
    's035';
    's036';
    's038';
    's041';
    %     %    's045';
    %     %    's046';
    's047';
    %     %    's048';
    's050';
    's051';
    's054';
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
for subj_idx=1:length(orig_subject)
    
    file_electrode=sprintf('electrode_tal_mri_%s_%s_061721.mat',roi_stem,orig_subject{subj_idx});
        
    load(sprintf('%s/%s/analysis/%s',root_path,orig_subject{subj_idx},file_electrode));
    electrode(subj_idx)=E;
    

end;

%load label
if(~isempty(file_label))
    [ll]=inverse_read_label(file_label);
else
    ll=[];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(~isempty(electrode))
    [orig_vertex_coords, orig_faces] = read_surf(sprintf('%s/%s/surf/%s.orig',subjects_dir,subject,hemi));
    [vertex_coords, faces] = read_surf(sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,hemi,surf));
    
    count=1;
    for e_idx=1:length(electrode)
        for c_idx=1:electrode(e_idx).n_contact
            
            aux2_point_coords(count,:)=electrode(e_idx).coord(c_idx,:);
            
            aux2_point_individual_color(count,:)=selected_electrode_color;
            if(c_idx==electrode(e_idx).contact_idx)
                aux2_point_individual_color(count,:)=selected_contact_color;
            end;
            
            aux2_point_name{count}=sprintf('%s_%d',electrode(e_idx).name, c_idx);;
            count=count+1;
        end;
    end;
end;




etc_render_fsbrain('overlay_cmap',[1.0000    0.4118    0.1608].*0.8,'overlay_smooth',[],'overlay_vertex',ll,'overlay_value',ones(length(ll),1),'overlay_threshold',[0.2 0.3],'curv_pos_color',[1 1 1].*0.5,'curv_neg_color',[1 1 1].*0.5,'alpha',0.3,'surf',surf,'hemi',hemi,'subject',subject,'vol',mri,'talxfm',(talxfm),'topo_aux2_point_coords',aux2_point_coords,'electrode',electrode);
view(90,30);
hgexport(gcf,sprintf('all_label_tal_brain_061721_%s',hemi), hgexport('factorystyle'),'Format','png');

global etc_render_fsbrain;
etc_render_fsbrain.aux2_point_individual_color=aux2_point_individual_color;
etc_render_fsbrain.selected_electrode_color=selected_electrode_color;
etc_render_fsbrain.selected_contact_color=selected_contact_color;
etc_render_fsbrain.selected_electrode_flag=0;
etc_render_fsbrain.selected_contact_flag=0;

etc_render_fsbrain_handle('redraw');


hgexport(gcf,sprintf('all_electrode_label_tal_brain_061721_%s',hemi), hgexport('factorystyle'),'Format','png');


return;
