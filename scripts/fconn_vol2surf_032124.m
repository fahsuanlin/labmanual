close all; clear all;



file_register='/Users/fhlin/workspace/seeg/s026/fmri_data/unpack/register_id.dat';
%provide the registration to structural freesurfer data


file_output={
'hippo_fconn_032124';
};

TR=2.0; %second

source_subject='s026';
target_subject='fsaverage';

file_in={
'/Users/fhlin/workspace/seeg/s026/fmri_analysis/hippo_fconn_native_vol_aseg_032124_gavg_hippo_left-anat.nii';
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdir=pwd;
for d_idx=1:length(file_in)

	for idx=1:length(file_output)
	%do this outside matlab....
	%make sure freesurfer environment, register file, and subjects directory are all set.

        [pathname,filename,ext]=fileparts(file_in{idx});
		fn0=sprintf('%s',file_in{idx});
		fn1=sprintf('%s/%s-lh.mgh',pathname,filename);
		fn2=sprintf('%s/%s-rh.mgh',pathname,filename);
		eval(sprintf('!mri_vol2surf --fwhm 10 --src %s --srcreg %s --hemi lh --noreshape --out %s',fn0,file_register,fn1));
		eval(sprintf('!mri_vol2surf --fwhm 10 --src %s --srcreg %s --hemi rh --noreshape --out %s',fn0,file_register,fn2));

		brain_lh = MRIread(fn1);
		fn3=sprintf('%s/%s-lh.stc',pathname,filename);
		stc=squeeze(brain_lh.vol); if(min(size(stc))==1) stc=stc'; end;
		inverse_write_stc(stc,[0:brain_lh.nvoxels-1],0,TR.*1e3,fn3);

		brain_rh = MRIread(fn2);
		fn4=sprintf('%s/%s-rh.stc',pathname,filename);
		stc=squeeze(brain_rh.vol); if(min(size(stc))==1) stc=stc'; end;
		inverse_write_stc(stc,[0:brain_rh.nvoxels-1],0,TR.*1e3,fn4);

		%morphing
        	fn_in=fn3;
            [stc,v,a,b]=inverse_read_stc(fn_in);
            inverse_write_stc(repmat(stc(:),[1 5]),v,a,b,'tmp-lh.stc');
        	fn_out=sprintf('%s_2_%s_%s',source_subject,target_subject,file_output{idx});
%        	cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'lh');
         	cmd=sprintf('!mne_make_movie --subject %s --stcin tmp --morph %s --stc %s --%s --smooth 5', source_subject, target_subject, fn_out, 'lh');
        	eval(cmd);

        	fn_in=fn4;        
            [stc,v,a,b]=inverse_read_stc(fn_in);
            inverse_write_stc(repmat(stc(:),[1 5]),v,a,b,'tmp-rh.stc');
		    fn_out=sprintf('%s_2_%s_%s',source_subject,target_subject,file_output{idx});
%        	cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'rh');
         	cmd=sprintf('!mne_make_movie --subject %s --stcin tmp --morph %s --stc %s --%s --smooth 5', source_subject, target_subject, fn_out, 'lh');
        	eval(cmd);


		eval(sprintf('!rm %s %s tmp-lh.stc tmp-rh.stc',fn1,fn2));
	end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pdir);
fprintf('DONE!\n');
