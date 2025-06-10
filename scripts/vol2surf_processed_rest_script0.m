close all; clear all;

subject={
'100206';
'100307';
'100408';
'100610';
'101006';
'101107';
'101309';
'101410';
'101915';
'102008';
%'102109';
%'102311';
%'102513';
%'102614';
%'102715';
%'102816';
%'103010';
%'103111';
%'103212';
%'103414';
%'103515';
%'103818';
%'104012';
%'104416';
%'104820';
%'105014';
%'105115';
%'105216';
%'105620';
%'105923';
%'106016';
%'106319';
%'106521';
%'106824';
%'107018';
%'107220';
%'107321';
%'107422';
%'107725';
%'108020';
};


file_stem={
'rfMRI_REST1_LR_hp2000_clean';
%'3T_rfMRI_REST1_LR';
%'3T_tfMRI_EMOTION_LR';
%'3T_tfMRI_SOCIAL_LR';
};

TR=0.72; %second

target_subject='fsaverage';

root_path='/space_lin1/hcp';

setenv('SUBJECTS_DIR','/space_lin1/hcp/subjects');



for subj_idx=1:length(subject)

	file_register=sprintf('/space_lin1/hcp/analysis/%s_register.dat',subject{subj_idx});

	if(~exist(sprintf('%s/%s/analysis',root_path,subject{subj_idx})))
		eval(sprintf('!mkdir %s/%s/analysis',root_path,subject{subj_idx}));
        end;

       for idx=1:length(file_stem)
        %do this outside matlab....
        %make sure freesurfer environment, register file, and subjects directory are all set.
                fn0=sprintf('%s/%s/Preprocessed/%s.nii.gz',root_path,subject{subj_idx},file_stem{idx});
		if(exist(fn0))
                fn1=sprintf('%s/%s/analysis/%s_%s-lh.mgh',root_path,subject{subj_idx},subject{subj_idx},file_stem{idx});
                fn2=sprintf('%s/%s/analysis/%s_%s-rh.mgh',root_path,subject{subj_idx},subject{subj_idx},file_stem{idx});
                eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi lh --noreshape --out %s',fn0,file_register,fn1));
                eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi rh --noreshape --out %s',fn0,file_register,fn2));

                brain_lh = MRIread(fn1);
                fn3=sprintf('%s/%s/analysis/%s_%s-lh.stc',root_path,subject{subj_idx},subject{subj_idx},file_stem{idx});
                stc=squeeze(brain_lh.vol); if(min(size(stc))==1) stc=stc'; end;
                inverse_write_stc(stc,[0:brain_lh.nvoxels-1],0,TR.*1e3,fn3);

                brain_rh = MRIread(fn2);
                fn4=sprintf('%s/%s/analysis/%s_%s-rh.stc',root_path,subject{subj_idx},subject{subj_idx},file_stem{idx});
                stc=squeeze(brain_rh.vol); if(min(size(stc))==1) stc=stc'; end;
                inverse_write_stc(stc,[0:brain_rh.nvoxels-1],0,TR.*1e3,fn4);

                %morphing
                fn_in=fn3;
                fn_out=sprintf('%s/%s/analysis/%s_2_%s_%s',root_path,subject{subj_idx},subject{subj_idx},target_subject,file_stem{idx});
                cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject{subj_idx}, fn_in, target_subject, fn_out, 'lh');
                eval(cmd);

                fn_in=fn4;
                fn_out=sprintf('%s/%s/analysis/%s_2_%s_%s',root_path,subject{subj_idx},subject{subj_idx},target_subject,file_stem{idx});
                cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject{subj_idx}, fn_in, target_subject, fn_out, 'rh');
                eval(cmd);


                eval(sprintf('!rm %s %s %s %s',fn1,fn2,fn3,fn4));
		end;
        end;

end;
