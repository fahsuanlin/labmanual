close all; clear all;

fmri_file={
'rfMRI_REST1_LR_hp2000_clean.nii.gz';
'rfMRI_REST1_RL_hp2000_clean.nii.gz';
'rfMRI_REST2_LR_hp2000_clean.nii.gz';
'rfMRI_REST2_RL_hp2000_clean.nii.gz';
};

setenv('SUBJECTS_DIR','/space_lin1/hcp/subjects');

subject='';
d=textread('subject_list_all.txt');
d=d(1:1);
for d_idx=1:length(d)
        subject{d_idx}=num2str(d(d_idx));
end;

for fmri_idx=1:length(fmri_file)

	[dummy, fmri_stem]=fileparts(fmri_file{fmri_idx});

	for subj_idx=1:length(subject)
%for subj_idx=1:length(subject)
		file_register=sprintf('/space_lin1/hcp/analysis/%s_%s_register.dat',subject{subj_idx},fmri_stem);

		if(~exist(file_register))
			fmri=sprintf('/space_lin1/hcp/%s/Preprocessed/%s',subject{subj_idx},fmri_file{fmri_idx});
			fmri_index=exist(fmri);
		
			if(sum(fmri_index)>eps)
				cmd=sprintf('!fslregister --s %s --mov  /space_lin1/hcp/%s/Preprocessed/%s --reg /space_lin1/hcp/analysis/%s_%s_register.dat --maxangle 70 --initxfm',subject{subj_idx},subject{subj_idx},fmri_file{fmri_idx},subject{subj_idx},fmri_stem);
				eval(cmd);
			end;

		else
			fprintf('subject [%s] has registration file already!\n',subject{subj_idx});
		end;
	end;
end;
