close all; clear all;

file_stem_stc={
    'seeg_wb_mne_091019_a_mne';
%    'seeg_wb_mne_091019_v_mne';
%    'seeg_wb_mne_091019_av_mne';
    };

subject='s031';
target_subject='fsaverage';


for f_idx=1:length(file_stem_stc)
    fstem=file_stem_stc{f_idx};
    fprintf('morphing [%s]...\n',fstem);
    for h=1:2
        if(h==1) hemi='lh'; else hemi='rh'; end;
        fn_in=sprintf('%s-%s.stc',fstem,hemi);
        fn_out=sprintf('%s_2_%s_%s',subject,target_subject,fstem);
        inverse_morph_stc(subject,fn_in,'hemi',hemi,'file_archive',fn_out);
        %cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_in, target_subject, fn_out, hemi);
        %eval(cmd);
    end;
end;