close all; clear all;

acc_mat={
'mb_run_1_acc.mat';
'mb_run_2_acc.mat';
};

acc_template={
'smsini_mb_run_1_ref.mgh'
'smsini_mb_run_2_ref.mgh'
};

file_register={
'register_01.dat';
'register_02.dat';
};

TR=2; %second

flag_morph=1;
flag_native=1;
flag_nii=1;

subject='060820';
target_subject='fsaverage';
output_stem='smsini';

for f_idx=1:length(acc_mat)
    load(acc_mat{f_idx});

    timeVec=[0:size(acc,4)-1].*TR;

    [dummy, fstem]=fileparts(acc_mat{f_idx});
    fstem=sprintf('%s_%s',output_stem, fstem);

    brain = MRIread(acc_template{f_idx});
    brain.vol=acc;
    brain.nframes=size(acc,4);
    MRIwrite(brain,sprintf('%s.mgh',fstem));

    if(flag_nii)
	        eval(sprintf('!mri_convert %s.mgh %s.nii',fstem,fstem));
    end;

    %do this outside matlab....
    %make sure freesurfer environment, register file, and subjects directory are all set.
    if(flag_morph)
	fn0=sprintf('%s.mgh',fstem);

	eval(sprintf('!mri_vol2vol --mov %s --reg %s --o %s-tal.2mm.mgh  --tal --talres 2',fn0,file_register{f_idx},fstem));
	if(flag_nii)
		fn1=sprintf('%s.nii',fstem);
		eval(sprintf('!mri_convert %s-tal.2mm.mgh %s-tal.2mm.nii',fstem,fstem));
	end;
    end;
    if(~flag_native)
	eval(sprintf('!rm %s',fn0));
	if(flag_nii)
		eval(sprintf('!rm %s',fn1));
	end;
    end;
end;



return;
