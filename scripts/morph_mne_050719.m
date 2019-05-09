close all; clear all;

file_stem_stc={
    'seeg_mne_050719_a_mne';
    'seeg_mne_050719_v_mne';
    'seeg_mne_050719_av_mne';
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
        cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_in, target_subject, fn_out, hemi);
        eval(cmd);
    end;
end;