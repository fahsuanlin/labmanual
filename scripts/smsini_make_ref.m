close all; clear all;

ref_mat={
'mb_run_1_ref.mat';
'mb_run_2_ref.mat';
'mb_run_3_ref.mat';
'mb_run_4_ref.mat';
'mb_run_5_ref.mat';
'mb_run_6_ref.mat';
};


for f_idx=1:length(ref_mat)
	load(ref_mat{f_idx});

	fmri_svbfile_fhlin(squeeze(ref),'gre_000.bfloat');
	[dummy,fstem]=fileparts(ref_mat{f_idx});
	%freesurfer convert
	cmd=sprintf('!mri_convert gre_000.bfloat %s.mgh -iid 1 0 0 -ijd 0 -1 0 -ikd 0 0 1 -iis 5 -ijs 5 -iks 5',fstem);
	eval(cmd);

	gre=MRIread(sprintf('%s.mgh',fstem));
	gre.vol=squeeze(gre.vol);
	MRIwrite(gre,sprintf('%s.mgh',fstem));

	eval('!rm gre_000.bfloat gre_000.hdr');

end;


return;


