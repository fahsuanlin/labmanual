close all; clear all;



subject='s063';

file_electrode_preop='electrode_042021_153243_s063';

electrode_select_name={
'x6';
't3';
};

electrode_select_color={
[0.6350 0.0780 0.1840];
[0.8500 0.3250 0.0980];
};

electrode_default_color=[.8 .9 .9];


flag_show_med=1;
flag_show_lat=1;
flag_show_electrode=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hemi={'lh'};


for hemi_idx=1:length(hemi)

    clear global etc_render_fsbrain;
    ff_med=figure;
    %set(ff_med,'visible','off');
    if(strcmp(hemi{hemi_idx},'lh'))
        etc_render_fsbrain('subject',subject,'hemi','lh','view_angle',[90,0],'alpha',0.4,'curv_neg_color',[1 1 1].*.6,'curv_pos_color',[1 1 1].*.6);
        view(90,0);
        view(-90,0);
    else
        etc_render_fsbrain('subject',subject,'hemi','rh','view_angle',[-90,0],'alpha',0.4,'curv_neg_color',[1 1 1].*.6,'curv_pos_color',[1 1 1].*.6);
        view(-90,0);
        view(90,0);
    end;

    if(flag_show_electrode)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the following show contact location
        global etc_render_fsbrain;


        load(file_electrode_preop);
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

                etc_render_fsbrain.aux2_point_individual_color(count,:)=electrode_default_color;

                for electrode_select_idx=1:length(electrode_select_name)
                    tmp=strcmp(lower(fn(2:end)),electrode_select_name{electrode_select_idx});
                    if(tmp)
                        etc_render_fsbrain.aux2_point_individual_color(count,:)=electrode_select_color{electrode_select_idx};
                    end;
                end;
                %etc_render_fsbrain.aux2_point_name{count}=sprintf('%s_%d',etc_render_fsbrain.electrode(e_idx).name, c_idx);;
                count=count+1;
            end;
        end;


        etc_render_fsbrain.selected_electrode_flag=0;
        etc_render_fsbrain.selected_contact_flag=0;

        etc_render_fsbrain_handle('redraw');

        %set(ff_med,'visible','off');
    end;


end;
