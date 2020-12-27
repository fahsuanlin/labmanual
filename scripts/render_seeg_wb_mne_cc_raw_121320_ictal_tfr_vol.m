close all; clear all;


file_stc={
    'seeg_wb_mne_cc_raw_121320_ictal_tfr15p0_dspm';
    };

subject='seeg027';

timeVec_range=[-100 2000];
hemi={'lh','rh'};

vol=MRIread('/space_lin2/fhlin/seeg_ictal_source/subjects/seeg027/mri/orig.mgz');

output_fstem=sprintf('seeg_wb_mne_cc_raw_121320_ictal_tfr15p0_dspm_vol');
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


talxfm=etc_read_xfm('file_xfm','/space_lin2/fhlin/seeg_ictal_source/subjects/seeg027/mri/transforms/talairach.xfm'); %for MAC/Linux

[vol_stc,vv,d0,d1,timeVec]=inverse_read_stc('seeg_wb_mne_cc_raw_121320_ictal_tfr15p0_dspm-vol.stc');

timeVec=timeVec(t_idx);
vol_stc=vol_stc(:,t_idx);

%sub-sampling in time
timeVec=timeVec(1:5:end);
vol_stc=vol_stc(:,1:5:end);

file_forward_mat='seeg_fwd_wb_091019.mat';
load(file_forward_mat);
for hemi_idx=1:2
%   file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hemi{hemi_idx},surf);
   file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hemi{hemi_idx},'orig');
   [A(hemi_idx).vertex_coords, A(hemi_idx).faces] = read_surf(file_surf);
end;

for time_idx=1:length(timeVec)
%    for time_idx=150:150
%     etc_render_fsbrain('subject',subject,'hemi','lh','overlay_value',squeeze(stc_lh(:,time_idx)),'overlay_vertex',v_lh,'overlay_threshold',[3 5],'overlay_smooth',5,'view_angle',[-90,0]);
%     hgexport(gcf,sprintf('img/%s_lh_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
%     im1=imread(sprintf('img/%s_rh_f%04d.png',output_fstem,time_idx));
%     close;
    
    etc_render_fsbrain('overlay_stc_timevec_idx',time_idx,'subject',subject,'hemi','rh','vol',vol,'talxfm',(talxfm),'overlay_vertex',v_rh,'overlay_vol_stc',vol_stc,'vol_A',A,'overlay_threshold',[8 12],'overlay_smooth',5,'view_angle',[-90,-20],'camposition_l',1500,'pt',[100 144 116]);
    figure(2);
    hgexport(gcf,sprintf('img/%s_vol_f%04d',output_fstem,time_idx),hgexport('factorystyle'),'Format','png');
    im2=imread(sprintf('img/%s_vol_f%04d.png',output_fstem,time_idx));
    close all;
    
    image(im2);
    axis off image;
    set(gcf,'pos',[680         315        800         800]);
    h=text(710, 750,sprintf('%0.0f ms',timeVec(time_idx)));
    set(h,'fontname','helvetica','fontsize',20,'color','c');
    hgexport(gcf,sprintf('img/tmp'),hgexport('factorystyle'),'Format','png');
    im=imread(sprintf('img/tmp.png'));
    im=imcrop(im,[120   171   600   440]);
    
    writeVideo(vidObj,im);
    close;
end;
close(vidObj);
