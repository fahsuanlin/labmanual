close all; clear all;


file_stc={
    'seeg_wb_mne_cc_fus_target_e4_tfr_raw_032822_e3_e4_stim3p0_freq006hz_dspm_sos';
    };

threshold={
    [6 10];
    };

subject='s052';

file_electrode_preop='electrode_042621_161326_s052';

flag_show_med=1;
flag_show_lat=1;
flag_show_electrode=1;

surface_coord=[1.50 -35.56  29.40]; %this is the surface coordinate to determine sagittal/coronal/axial slices
target_timeVec=[-50:25:500];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baseline=[-200 -20];

%timeVec_range=[-200 800];

hemi={'lh'};

trial=[1];

load('erf_raw_fus_target_e4_032322.mat');

hh={'lh','rh'};
for hemi_idx=1:2
   file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hh{hemi_idx},'orig');
   [A(hemi_idx).vertex_coords, A(hemi_idx).faces] = read_surf(file_surf);
end;


for stc_idx=1:length(file_stc)
    for trial_idx=1:length(trial)
         [stc_lh,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s-lh.stc',file_stc{stc_idx}));
         [stc_rh,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s-rh.stc',file_stc{stc_idx}));
        

        im2_med=[];
        im2_lat=[];
        
        
        output_fstem=sprintf('%s_trial%03d_vol',file_stc{stc_idx},trial(trial_idx));
        
        for hemi_idx=1:length(hemi)
            
            clear global etc_render_fsbrain;
            ff_med=figure; set(ff_med,'visible','off');
            if(strcmp(hemi{hemi_idx},'lh'))
                etc_render_fsbrain('subject',subject,'hemi','lh','overlay_stc',stc_lh,'overlay_vertex',v_lh,'vol_A',A,'overlay_threshold',threshold{stc_idx},'overlay_smooth',5,'view_angle',[90,0]);
                view(90,0);
                view(-90,0);
                etc_render_fsbrain_handle('draw_pointer','surface_coord',surface_coord);

            else
                etc_render_fsbrain('subject',subject,'hemi','rh','overlay_stc',stc_rh,'overlay_vertex',v_rh,'vol_A',A,'overlay_threshold',threshold{stc_idx},'overlay_smooth',5,'view_angle',[-90,0]);
                view(-90,0);
                view(90,0);
                etc_render_fsbrain_handle('draw_pointer','surface_coord',surface_coord);
                
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
                        
                        %etc_render_fsbrain.aux2_point_name{count}=sprintf('%s_%d',etc_render_fsbrain.electrode(e_idx).name, c_idx);;
                        count=count+1;
                    end;
                end;
                
                etc_render_fsbrain_handle('redraw');
                set(ff_med,'visible','off');
            end;
            

            [~, timeVec_idx] = min( abs(timeVec(:) - target_timeVec(:).'), [], 1 );
            for time_idx=1:length(timeVec_idx)
                fprintf('rendering <%s>:: [%0d|%04d] (timeVec=%1.1f; %2.2f%%)...\r',output_fstem, time_idx,length(timeVec_idx),timeVec(timeVec_idx(time_idx)),time_idx./length(timeVec_idx).*100);

                etc_render_fsbrain_update_time(timeVec_idx(time_idx));
                etc_render_fsbrain_handle('redraw');
                etc_render_fsbrain_handle('draw_pointer','surface_coord',surface_coord);
                figure(etc_render_fsbrain.fig_vol);
                hgexport(gcf,sprintf('img/%s_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
            end;
        end;        

        if(ismac)
            vidObj = VideoWriter(output_fstem,'MPEG-4');
        else
            vidObj = VideoWriter(output_fstem,'Archival');
        end;
        vidObj.FrameRate=24; %fps
        open(vidObj);
        
        for time_idx=1:length(timeVec_idx)
            for hemi_idx=1:length(hemi)
                fprintf('create animation frame <%s>::<<%s>> [%0d|%04d] (timeVec=%1.1f; %2.2f%%)...\r',output_fstem, hemi{hemi_idx},time_idx,length(timeVec_idx),timeVec(timeVec_idx(time_idx)),time_idx./length(timeVec_idx).*100);
                
                tmp=imread(sprintf('img/%s_f%04d.png',output_fstem,time_idx));
            end;
            
            ff_img=figure;
            image(tmp);
            axis off image;
            set(gca,'pos',[0 0 1 1]);
            set(gcf,'pos',[200         800        size(tmp,2)*2         size(tmp,1)*2]);
            h=text(650, 460, sprintf('%0.0f ms',timeVec(timeVec_idx(time_idx))));
            set(h,'fontname','helvetica','fontsize',40,'color','b');
            hgexport(ff_img,sprintf('img/%s_tmp',output_fstem),hgexport('factorystyle'),'Format','png');
            im=imread(sprintf('img/%s_tmp.png',output_fstem));
            
            writeVideo(vidObj,im);
            close;
        end;        
        fprintf('\n');
        close(vidObj);
    end;
end;
