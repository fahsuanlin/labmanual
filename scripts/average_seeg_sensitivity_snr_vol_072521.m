close all; clear all;

subject={
    's026';
    's027';
    's031';
    's032';
    's033';
    's034';
    's036';
    's041';
    's045';
    's047';
    };


file_snr_stem={
	'seeg_wb_sensitivity_072521_snr-vol-tal.2mm.mgh';
        'seeg_wb_sensitivity_072521_s-vol-tal.2mm.mgh';
};
output_mne_stem={
    'average_seeg_wb_sensitiivty_072521_snr-vol-tal.2mm.mgh';
    'average_seeg_wb_sensitiivty_072521_s-vol-tal.2mm.mgh';
};

path_root='/space_lin2/fhlin/seeg';
%path_root='/Users/fhlin_admin/workspace/seeg';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare brain mask
%mri_overlay=MRIread('/Applications/freesurfer/subjects/fsaverage/mri/brainmask.mgz');
d=sprintf('%s/subjects/fsaverage/mri/brainmask.mgz',getenv('FREESURFER_HOME'));
mri_overlay=MRIread(d);

mask_subject='fsaverage';

%targ=MRIread('/Applications/freesurfer/average/mni305.cor.subfov2.mgz'); %MNI-Talairach space with 2mm resolution (for MAC)
d=sprintf('%s/average/mni305.cor.subfov2.mgz',getenv('FREESURFER_HOME'));
targ=MRIread(d);
%targ=MRIread(sprintf('%s/average/mni305.cor.subfov2.mgz',getenv('FREESURFER_HOME'))); %MNI-Talairach space with 2mm resolution (for server)

fprintf('loading transformation for subject p%s]...\n',mask_subject);
mov_xfm=etc_read_xfm('subject',mask_subject);

R=mri_overlay.tkrvox2ras*inv(mri_overlay.vox2ras)*inv(mov_xfm)*(targ.vox2ras)*inv(targ.tkrvox2ras);
mri_mask_tal=etc_MRIvol2vol(mri_overlay,targ,R);
dd=zeros(size(mri_mask_tal.vol));
dd(find(mri_mask_tal.vol(:)))=1;
mri_mask_tal.vol=dd;
MRIwrite(mri_mask_tal,sprintf('brainmask-tal.2mm.mgh'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for f_idx=1:length(file_snr_stem)
    for subj_idx=1:length(subject)
        
        fn=sprintf('%s/%s/analysis/%s',path_root,subject{subj_idx},file_snr_stem{f_idx});
        fprintf('loading [%s]...\n',fn);
        d=MRIread(fn);

        tmp_snr=d.vol;
        

        fprintf('smoothing');
        [dd,k]=fmri_smooth(tmp_snr(:,:,:,1),6,'vox',[2 2 2]);
        tmp_snr(:,:,:,1)=dd.*(mri_mask_tal.vol); %apply a mask
        snr_all(:,:,:,subj_idx)=tmp_snr(:,:,:,1);

    end;
    vol=mean(snr_all,4);

    
    d.vol=vol;
    ff=sprintf('%s',output_mne_stem{f_idx});
    fprintf('saving [%s]...\n',ff);
    MRIwrite(d,ff);
end;

return;
