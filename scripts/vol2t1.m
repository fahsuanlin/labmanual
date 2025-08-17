close all; clear all;

file_register={
'../../../../motor_analysis/register_006.dat';
'../../../../motor_analysis/register_009.dat';
};
%provide the registration to structural freesurfer data


file_output={
't1_sfmcprstc';
};

TR=2.0; %second

source_subject='s021';
target_subject='s021';
dirs={
'bold/006';
'bold/009';
};


%provide the output file name for InI reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdir=pwd;
for d_idx=1:length(dirs)
	dd=sprintf('%s/%s',pdir,dirs{d_idx});
	cd(dd);

	for idx=1:length(file_output)
	%do this outside matlab....
	%make sure freesurfer environment, register file, and subjects directory are all set.
		fn0=sprintf('%s.nii',file_output{idx});
		eval(sprintf('!mri_vol2vol --mov sfmcprstc.nii --reg %s --o %s.nii  --no-resample  --fstarg',file_register{d_idx},file_output{1}));

	end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pdir);
fprintf('DONE!\n');

