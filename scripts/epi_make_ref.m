close all; clear all;

ref_mat={
'../fmri_data/unpack/bold/005/f.nii';
'../fmri_data/unpack/bold/018/f.nii';
'../fmri_data/unpack/bold/027/f.nii';
'../fmri_data/unpack/bold/036/f.nii';
'../fmri_data/unpack/bold/045/f.nii';
'../fmri_data/unpack/bold/054/f.nii';
'../fmri_data/unpack/bold/063/f.nii';
};

output_stem='epi';

for f_idx=1:length(ref_mat)
	%load(ref_mat{f_idx});
    ref=MRIread(ref_mat{f_idx});
    ref=mean(ref.vol,4);
    for s_idx=1:size(ref,3)
        ref(:,:,s_idx)=rot90(ref(:,:,s_idx));
    end;

	fmri_svbfile_fhlin(squeeze(ref),'gre_000.bfloat');
	[dummy,fstem]=fileparts(ref_mat{f_idx});

    run_stem=regexp(dummy,'\-?[\d\.]+$','match','once');
	fstem=sprintf('%s_%s_%s',output_stem, run_stem, fstem);
	%freesurfer convert
	cmd=sprintf('!mri_convert gre_000.bfloat %s.mgh -iid 1 0 0 -ijd 0 -1 0 -ikd 0 0 1 -iis 3.28125 -ijs 3.28125 -iks 5 ',fstem);
	eval(cmd);

	gre=MRIread(sprintf('%s.mgh',fstem));
	gre.vol=squeeze(gre.vol);
	MRIwrite(gre,sprintf('%s.mgh',fstem));

	eval('!rm gre_000.bfloat gre_000.hdr');

end;


return;
