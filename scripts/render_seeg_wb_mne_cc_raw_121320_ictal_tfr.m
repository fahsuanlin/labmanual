close all; clear all;


file_stc={
    'seeg_wb_mne_cc_raw_121320_ictal_tfr15p0_dspm';
    };

subject='seeg027';

timeVec_range=[-100 2000];
hemi={'lh','rh'};


output_fstem=sprintf('seeg_wb_mne_cc_raw_121320_ictal_tfr15p0_dspm');
vidObj = VideoWriter(output_fstem,'MPEG-4');
vidObj.FrameRate=16; %fps
open(vidObj);

[stc_lh,v_lh,a,b,timeVec]=inverse_read_stc(sprintf('%s-lh.stc',file_stc{1}));
[stc_rh,v_rh,a,b,timeVec]=inverse_read_stc(sprintf('%s-rh.stc',file_stc{1}));

t_idx=find((timeVec>=min(timeVec_range))&(timeVec<=max(timeVec_range)));
timeVec=timeVec(t_idx);
stc_lh=stc_lh(:,t_idx);
stc_rh=stc_rh(:,t_idx);

im1=[];
im2=[];
for time_idx=1:length(timeVec)
%     etc_render_fsbrain('subject',subject,'hemi','lh','overlay_value',squeeze(stc_lh(:,time_idx)),'overlay_vertex',v_lh,'overlay_threshold',[3 5],'overlay_smooth',5,'view_angle',[-90,0]);
%     hgexport(gcf,sprintf('img/%s_lh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
%     im1=imread(sprintf('img/%s_rh_f%04d.png',output_fstem,time_idx));
%     close;
    
    etc_render_fsbrain('subject',subject,'hemi','rh','overlay_value',squeeze(stc_rh(:,time_idx)),'overlay_vertex',v_rh,'overlay_threshold',[8 12],'overlay_smooth',5,'view_angle',[-90,-20],'camposition_l',1500);
    hgexport(gcf,sprintf('img/%s_rh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
    im2=imread(sprintf('img/%s_rh_f%04d.png',output_fstem,time_idx));
    close;
    
    image(cat(2,im1,im2));
    axis off image;
    set(gcf,'pos',[680         315        800         800]);
    h=text(40, 340,sprintf('%0.0f ms',timeVec(time_idx)));
    set(h,'fontname','helvetica','fontsize',20,'color','b');
    hgexport(gcf,sprintf('img/tmp'),hgexport('factorystyle'),'Format','png');
    im=imread(sprintf('img/tmp.png'));
    
    writeVideo(vidObj,im);
    close;
end;
close(vidObj);
