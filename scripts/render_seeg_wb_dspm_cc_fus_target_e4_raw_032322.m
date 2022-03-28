close all; clear all;


file_stc={
    'seeg_wb_mne_cc_fus_target_e4_raw_032322_fus_stim_mne_sos_dspm';
    'seeg_wb_mne_cc_fus_target_e4_raw_032322_post_stim_mne_sos_dspm';
    };

threshold={
    [3 10];
    [3 10];
    };

subject='s052';

file_electrode_preop='electrode_042621_161326_s052';

electrode_select_name={
    'e4';
    'f9';
    'g1';
    };

electrode_select_color={
    [0.6350 0.0780 0.1840];
    [0.8500 0.3250 0.0980];
    [0.4940 0.1840 0.5560];
    };

electrode_default_color=[.8 .9 .9];


flag_show_med=1;
flag_show_lat=1;
flag_show_electrode=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baseline=[-200 -20];

%timeVec_range=[-200 800];

hemi={'lh'};

trial=[1];

load('erf_raw_fus_target_e4_032322.mat');

for stc_idx=1:length(file_stc)
    for trial_idx=1:length(trial)
%         [stc_lh_x,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_x-lh.stc',file_stc{stc_idx},trial(trial_idx)));
%         [stc_rh_x,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_x-rh.stc',file_stc{stc_idx},trial(trial_idx)));
%         [stc_lh_y,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_y-lh.stc',file_stc{stc_idx},trial(trial_idx)));
%         [stc_rh_y,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_y-rh.stc',file_stc{stc_idx},trial(trial_idx)));
%         [stc_lh_z,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_z-lh.stc',file_stc{stc_idx},trial(trial_idx)));
%         [stc_rh_z,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s_trial%03d_z-rh.stc',file_stc{stc_idx},trial(trial_idx)));
%         stc_lh=sqrt(stc_lh_x.^2+stc_lh_y.^2+stc_lh_z.^2);
%         stc_rh=sqrt(stc_rh_x.^2+stc_rh_y.^2+stc_rh_z.^2);
         [stc_lh,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s-lh.stc',file_stc{stc_idx}));
         [stc_rh,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s-rh.stc',file_stc{stc_idx}));
        

%         t_idx=find((timeVec>=min(timeVec_range))&(timeVec<=max(timeVec_range)));
%         timeVec=timeVec(t_idx);
%         stc_lh=stc_lh(:,t_idx);
%         stc_rh=stc_rh(:,t_idx);

        im2_med=[];
        im2_lat=[];
        
%         timeVec=timeVec(1:20:end);
%         stc_lh=stc_lh(:,1:20:end);
%         stc_rh=stc_rh(:,1:20:end);
        
        output_fstem=sprintf('%s_trial%03d',file_stc{stc_idx},trial(trial_idx));
        
        for hemi_idx=1:length(hemi)
            
            clear global etc_render_fsbrain;
            ff_med=figure; set(ff_med,'visible','off');
            if(strcmp(hemi{hemi_idx},'lh'))
                etc_render_fsbrain('subject',subject,'hemi','lh','overlay_stc',stc_lh,'overlay_vertex',v_lh,'overlay_threshold',threshold{stc_idx},'overlay_smooth',5,'view_angle',[90,0],'alpha',0.4,'curv_neg_color',[1 1 1].*.6,'curv_pos_color',[1 1 1].*.6);
                view(90,0);
                view(-90,0);
            else
                etc_render_fsbrain('subject',subject,'hemi','rh','overlay_stc',stc_rh,'overlay_vertex',v_rh,'overlay_threshold',threshold{stc_idx},'overlay_smooth',5,'view_angle',[-90,0],'alpha',0.4,'curv_neg_color',[1 1 1].*.6,'curv_pos_color',[1 1 1].*.6);
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
                
                set(ff_med,'visible','off');
            end;
            
            for time_idx=1:length(timeVec)
                fprintf('rendering <%s>:: [%0d|%04d] (timeVec=%1.1f; %2.2f%%)...\r',output_fstem, time_idx,length(timeVec),timeVec(time_idx),time_idx./length(timeVec).*100);
                
                if(flag_show_med)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % medial view
                    etc_render_fsbrain_update_time(time_idx);
                    if(strcmp(hemi{hemi_idx},'lh'))
                        %view(90,0);
                        view(109, 36);
                        %hgexport(gcf,sprintf('img/%s_med_lh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
                        print(sprintf('img/%s_med_lh_f%04d',output_fstem,time_idx),'-dpng','-r300');      
                    else
                        view(-90,0);                        
                        %hgexport(gcf,sprintf('img/%s_med_rh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
                        print(sprintf('img/%s_med_rh_f%04d',output_fstem,time_idx),'-dpng','-r300');   
                    end;
                end;
                
                if(flag_show_lat)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % lateral view
                    etc_render_fsbrain_update_time(time_idx);
                    if(strcmp(hemi{hemi_idx},'lh'))
                        %view(-90,0);
                        view(-109,36);
                        %hgexport(gcf,sprintf('img/%s_lat_lh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
                        print(sprintf('img/%s_lat_lh_f%04d',output_fstem,time_idx),'-dpng','-r300');      
                    else
                        view(90,0);                        
                        %hgexport(gcf,sprintf('img/%s_lat_rh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
                        print(sprintf('img/%s_lat_rh_f%04d',output_fstem,time_idx),'-dpng','-r300');      
                    end;
                end;
            end;
        end;        

        if(ismac)
            vidObj = VideoWriter(output_fstem,'MPEG-4');
        else
            vidObj = VideoWriter(output_fstem,'Archival');
        end;
        vidObj.FrameRate=24; %fps
        open(vidObj);
        
        for time_idx=1:length(timeVec)
            for hemi_idx=1:length(hemi)
                fprintf('create animation frame <%s>::<<%s>> [%0d|%04d] (timeVec=%1.1f; %2.2f%%)...\r',output_fstem, hemi{hemi_idx},time_idx,length(timeVec),timeVec(time_idx),time_idx./length(timeVec).*100);
                im2_lh_lat=[];
                im2_rh_lat=[];
                im2_lh_med=[];
                im2_rh_med=[];
                
                if(strcmp(hemi{hemi_idx},'lh'))
                    if(flag_show_lat)
                        im2_lh_lat=imread(sprintf('img/%s_lat_lh_f%04d.png',output_fstem,time_idx));
                    end;
                    if(flag_show_med)
                        im2_lh_med=imread(sprintf('img/%s_med_lh_f%04d.png',output_fstem,time_idx));
                    end;
                end;
                if(strcmp(hemi{hemi_idx},'rh'))
                    if(flag_show_lat)
                        im2_rh_lat=imread(sprintf('img/%s_lat_rh_f%04d.png',output_fstem,time_idx));
                    end;
                    if(flag_show_med)
                        im2_rh_med=imread(sprintf('img/%s_med_rh_f%04d.png',output_fstem,time_idx));
                    end;
                end;
            end;
            
            ff_img=figure; set(ff_img,'visible','off');
            tmp=cat(2,im2_lh_lat,im2_lh_med,im2_rh_med,im2_rh_lat);
            image(tmp);
            axis off image;
            set(gca,'pos',[0 0 1 1]);
            h=text(120, 1400, sprintf('%1.0f Hz',erf_all(1).timeVec(time_idx)));
            set(h,'fontname','helvetica','fontsize',20,'color','b');
            %hgexport(ff_img,sprintf('img/%s_tmp',output_fstem),hgexport('factorystyle'),'Format','png');
            print(sprintf('img/%s_tmp.png',output_fstem),'-dpng','-r300')
            im=imread(sprintf('img/%s_tmp.png',output_fstem));

            writeVideo(vidObj,im);
            close;
        end;        
        fprintf('\n');
        close(vidObj);
    end;
end;
