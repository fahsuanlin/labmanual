close all; clear all;


file_register='../../register.dat';
%provide the registration to structural freesurfer data


file_output={
'sfmcprstc';
};

TR=2.0; %second

source_subject='s026';
target_subject='fsaverage';
dirs={
%'bold/030';
'bold/032';
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
		fn1=sprintf('%s-lh.mgh',file_output{idx});
		fn2=sprintf('%s-rh.mgh',file_output{idx});
		eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi lh --noreshape --out %s',fn0,file_register,fn1));
		eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi rh --noreshape --out %s',fn0,file_register,fn2));

		brain_lh = MRIread(fn1);
		fn3=sprintf('%s-lh.stc',file_output{idx});
		stc=squeeze(brain_lh.vol); if(min(size(stc))==1) stc=stc'; end;
		inverse_write_stc(stc,[0:brain_lh.nvoxels-1],0,TR.*1e3,fn3);

		brain_rh = MRIread(fn2);
		fn4=sprintf('%s-rh.stc',file_output{idx});
		stc=squeeze(brain_rh.vol); if(min(size(stc))==1) stc=stc'; end;
		inverse_write_stc(stc,[0:brain_rh.nvoxels-1],0,TR.*1e3,fn4);

		%morphing
        	fn_in=fn3;
        	fn_out=sprintf('%s_2_%s_%s',source_subject,target_subject,file_output{idx});
        	cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'lh');
        	eval(cmd);

        	fn_in=fn4;        
		fn_out=sprintf('%s_2_%s_%s',source_subject,target_subject,file_output{idx});
        	cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'rh');
        	eval(cmd);


	%	eval(sprintf('!rm %s %s',fn1,fn2));
	end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pdir);
fprintf('DONE!\n');
