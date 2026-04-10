close all; clear all;


file_register={
    'bb_register_01.dat';
    };
%provide the registration to structural freesurfer data
file_fmri_stem='fmcprstc_topup';

file_output={
    'fmcprstc_topup_040';
    };

TR=2.5; %second

source_subject='hjlee';
target_subject='fsaverage';
dirs={
    '../fmri_data/unpack/bold/040';
    };


surf_name={
    'white';
    'gray01';
    'gray02';
    'gray03';
    'gray04';
    'gray05';
    'gray06';
    'gray07';
    'gray08';
    'gray09';
    'pial';
    };

volumeSmooth = 0;
surfaceSmooth = 5;

%provide the output file name for InI reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdir=pwd;
for d_idx=1:length(dirs)
    
    for surf_idx=1:length(surf_name)
        
            %do this outside matlab....
            %make sure freesurfer environment, register file, and subjects directory are all set.
            fn0=sprintf('%s/%s.nii.gz',dirs{d_idx},file_fmri_stem);
            fn1=sprintf('%s_%s-lh.mgh',file_fmri_stem,surf_name{surf_idx});
            fn2=sprintf('%s_%s-rh.mgh',file_fmri_stem,surf_name{surf_idx});
            %eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi lh --noreshape --out %s',fn0,file_register,fn1));
            %eval(sprintf('!mri_vol2surf --icoorder 5 --fwhm 10 --src %s --srcreg %s --hemi rh --noreshape --out %s',fn0,file_register,fn2));
            
            eval(sprintf('!mri_vol2surf --interp nearest --fwhm %f --surf-fwhm %f --src %s --srcreg %s --hemi lh --surf %s --out %s',volumeSmooth,surfaceSmooth,fn0,file_register{d_idx},surf_name{surf_idx},fn1));
            eval(sprintf('!mri_vol2surf --interp nearest --fwhm %f --surf-fwhm %f --src %s --srcreg %s --hemi rh --surf %s --out %s',volumeSmooth,surfaceSmooth,fn0,file_register{d_idx},surf_name{surf_idx},fn2));
            
            brain_lh = MRIread(fn1);
            fn3=sprintf('%s_%s-lh.stc',file_output{d_idx},surf_name{surf_idx});
            stc=squeeze(brain_lh.vol); if(min(size(stc))==1) stc=stc'; end;
            inverse_write_stc(stc,[0:brain_lh.nvoxels-1],0,TR.*1e3,fn3);
            
            brain_rh = MRIread(fn2);
            fn4=sprintf('%s_%s-rh.stc',file_output{d_idx},surf_name{surf_idx});
            stc=squeeze(brain_rh.vol); if(min(size(stc))==1) stc=stc'; end;
            inverse_write_stc(stc,[0:brain_rh.nvoxels-1],0,TR.*1e3,fn4);
            
            %morphing
            fn_in=fn3;
            fn_out=sprintf('%s_2_%s_%s_%s',source_subject,target_subject,file_output{d_idx},surf_name{surf_idx});
            cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'lh');
            eval(cmd);
            
            fn_in=fn4;
            fn_out=sprintf('%s_2_%s_%s_%s',source_subject,target_subject,file_output{d_idx},surf_name{surf_idx});
            cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', source_subject, fn_in, target_subject, fn_out, 'rh');
            eval(cmd);
            
            
            %eval(sprintf('!rm %s %s %s %s',fn1,fn2,fn3,fn4));
            eval(sprintf('!rm %s %s %s %s',fn1,fn2));
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('DONE!\n');
