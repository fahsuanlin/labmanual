close all; clear all;

tsnr_stc={
    'tsnr_20ch_tsnr';
    'tsnr_32ch_tsnr';
};

threshold=[30 100];


for f_idx=1:length(tsnr_stc)

    [dummy, output_stem]=fileparts(tsnr_stc{f_idx});
    
    h=etc_render_fsbrain_stc({tsnr_stc{f_idx}},threshold,'overlay_smooth',5);
    %etc_render_fsbrain_stc({fstem{f_idx}},[Zid_l Zid_h],'overlay_smooth',10,'overlay_exclude_fstem','exclude','surf','pial','flag_curv',0, 'default_solid_color',[1 1 1].*.8);
    set(gcf,'pos',[680   770   960   400]);
    hgexport(gcf,output_stem, hgexport('factorystyle'),'Format','png');
    
end;