close all; clear all;

file_stc={
    'average_041620_seeg_wb_mne_091019_l_z';
    'average_041620_seeg_wb_mne_091019_n_z';
    };

output_fstem={
    'average_041620_seeg_wb_mne_091019_l_z';
    'average_041620_seeg_wb_mne_091019_n_z';
    };

threshold={
    [30 70];
    [30 70];
    };

for f_idx=1:length(file_stc)
    [stc,v,a,b,timeVec]=inverse_read_stc(sprintf('%s-lh.stc',file_stc{f_idx}));
    timeVec_idx=[1:length(timeVec)];
    timeVec_idx=timeVec_idx(1:20:end); %about every 10 ms

    
    fstem=output_fstem{f_idx};
    
    vidObj = VideoWriter(output_fstem{f_idx},'MPEG-4');
    vidObj.FrameRate=8; %fps
    open(vidObj);
    
    for t_idx=1:length(timeVec_idx)
        etc_render_fsbrain_stc({file_stc{f_idx}},threshold{f_idx},'time_idx',timeVec_idx(t_idx),'flag_view_default4',1);        
        
        hgexport(gcf,sprintf('img/%s_%04d',fstem,t_idx),hgexport('factorystyle'),'Format','png');
        close;
        
        im=imread(sprintf('img/%s_%04d.png',fstem,t_idx));
        %im=imcrop(im,[226   174   168   168]);
        image(im);
        axis off image;
        set(gcf,'pos',[1275         801        1800         600]);
        h=text(480,222,sprintf('%1.0f ms',timeVec(timeVec_idx(t_idx))));
        set(h,'fontname','helvetica','fontsize',20,'color','k');
        hgexport(gcf,sprintf('img/tmp'),hgexport('factorystyle'),'Format','png');
        im=imread(sprintf('img/tmp.png'));
        
        writeVideo(vidObj,im);
        close;
    end;
    close(vidObj);
end;
