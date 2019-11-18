close all; clear all;

acc_mat={
'mb_run_1_acc.mat';
'mb_run_2_acc.mat';
'mb_run_3_acc.mat';
'mb_run_4_acc.mat';
'mb_run_5_acc.mat';
'mb_run_6_acc.mat';
};

acc_template={
'mb_run_1_ref.mgh';
'mb_run_2_ref.mgh';
'mb_run_3_ref.mgh';
'mb_run_4_ref.mgh';
'mb_run_5_ref.mgh';
'mb_run_6_ref.mgh';
};

file_register={
'register_01.dat';
'register_02.dat';
'register_03.dat';
'register_04.dat';
'register_05.dat';
'register_06.dat';
};

TR=0.1; %second

flag_morph=1;
subject='191115_social';
target_subject='fsaverage';

for f_idx=1:length(acc_mat)
        load(acc_mat{f_idx});

	timeVec=[0:size(acc,4)-1].*TR;

	[dummy, fstem]=fileparts(acc_mat{f_idx});

    brain = MRIread(acc_template{f_idx});
    brain.vol=acc;
    brain.nframes=size(acc,4);
    MRIwrite(brain,sprintf('%s.mgh',fstem));

    %do this outside matlab....
    %make sure freesurfer environment, register file, and subjects directory are all set.
    fn0=sprintf('%s.mgh',fstem);
    fn1=sprintf('%s-lh.mgh',fstem);
    fn2=sprintf('%s-rh.mgh',fstem);
    eval(sprintf('!mri_vol2surf --fwhm 10 --src %s --srcreg %s --hemi lh --noreshape --out %s',fn0,file_register{f_idx},fn1));
    eval(sprintf('!mri_vol2surf --fwhm 10 --src %s --srcreg %s --hemi rh --noreshape --out %s',fn0,file_register{f_idx},fn2));

    brain_lh = MRIread(fn1);
    fn3=sprintf('%s-lh.stc',fstem);
    stc=squeeze(brain_lh.vol);
    inverse_write_stc(stc,[0:brain_lh.nvoxels-1],timeVec(1).*1e3,mean(diff(timeVec)).*1e3,fn3);

    brain_rh = MRIread(fn2);
    fn4=sprintf('%s-rh.stc',fstem);
    stc=squeeze(brain_rh.vol);
    inverse_write_stc(stc,[0:brain_rh.nvoxels-1],timeVec(1).*1e3,mean(diff(timeVec)).*1e3,fn4);

    eval(sprintf('!rm %s %s %s',fn0,fn1,fn2));

    if(flag_morph)
        fprintf('morphing....\n');
        for h=1:2
            if(h==1) hemi='lh'; else hemi='rh'; end;
            fn_in=sprintf('%s-%s.stc',fstem,hemi);
            fn_out=sprintf('%s_2_%s_%s',subject,target_subject,fstem);
            cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_in, target_subject, fn_out, hemi);
            eval(cmd);
        end;
        %eval(sprintf('!rm %s %s',fn3,fn4));
    end;
end;



return;



