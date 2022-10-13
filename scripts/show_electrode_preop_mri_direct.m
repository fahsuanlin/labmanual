close all; clear all;

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/music_layer_seeg/subjects/'); %for MAC/Linux
%setenv('SUBJECTS_DIR','D:\fhlin\Users\fhlin\workspace\seeg\subjects'); %for PC

mri=MRIread('/Users/fhlin/workspace/music_layer_seeg/subjects/s035/mri/orig.mgz'); %for MAC/Linux
%mri=etc_MRIread('D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\orig.mgz'); %for PC

% do registration between pre- and post-OP by the following command:
%
% cd /Users/fhlin/workspace/music_layer_seeg/subjects/s035_post/tmp
% bbregister --s s035 --mov ../mri/orig.mgz --init-fsl --reg register.dat --t1
%
%load the Talairach transformation matrix from the "pre-OP" data
talxfm=etc_read_xfm('file_xfm','/Users/fhlin/workspace/music_layer_seeg/subjects/s035/mri/transforms/talairach.xfm'); %for MAC/Linux
%talxfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036\mri\transforms\talairach.xfm'); %for PC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show brain


etc_render_fsbrain('surf','orig','hemi','rh','subject','s035','vol',mri,'talxfm',(talxfm),'alpha',0.5);
view(90,30);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show electrodes

file_mat='electrode_041819_114159_s035.mat'; %electrode coordinates in post-op MRI
load(file_mat);

electrode_select_color={
    [0.6350 0.0780 0.1840];
    };

global etc_render_fsbrain;

%update electrode contact coordinates
etc_render_fsbrain.electrode=electrode;

etc_render_fsbrain.aux2_point_coords=[];
etc_render_fsbrain.aux2_point_name={};
count=1;
for e_idx=1:length(etc_render_fsbrain.electrode)
    for c_idx=1:etc_render_fsbrain.electrode(e_idx).n_contact
        
        etc_render_fsbrain.aux2_point_coords(count,:)=etc_render_fsbrain.electrode(e_idx).coord(c_idx,:);
        
        if(strcmp(etc_render_fsbrain.surf,'orig')|strcmp(etc_render_fsbrain.surf,'smoothwm')|strcmp(etc_render_fsbrain.surf,'pial'))
            
        else
            %fprintf('surface <%s> not "orig"/"smoothwm"/"pial". Electrode contacts locations are updated to the nearest location of this surface.\n',etc_render_fsbrain.surf);
            
            tmp=etc_render_fsbrain.aux2_point_coords(count,:);
            
            vv=etc_render_fsbrain.orig_vertex_coords;
            dist=sqrt(sum((vv-repmat([tmp(1),tmp(2),tmp(3)],[size(vv,1),1])).^2,2));
            [min_dist,min_dist_idx]=min(dist);
            if(~isnan(min_dist))
                etc_render_fsbrain.aux2_point_coords(count,:)=etc_render_fsbrain.vertex_coords(min_dist_idx,:);
            end;
        end;
        
        fn=sprintf('%s%d',etc_render_fsbrain.electrode(e_idx).name,c_idx);
        
        etc_render_fsbrain.aux2_point_individual_color(count,:)=electrode_select_color{1};
        
        count=count+1;
    end;
end;


etc_render_fsbrain.selected_electrode_flag=1;
etc_render_fsbrain.selected_contact_flag=1;

etc_render_fsbrain_handle('redraw');
