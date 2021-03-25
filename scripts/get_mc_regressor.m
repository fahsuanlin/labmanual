close all; clear all;

fmc_stem='fmcpr.mcdat';

bold_dir='/space_lin2/fhlin/7t_music_skku/LAM_AUD_BHC_simple/unpack';

dirs={
'bold/007';
'bold/008';
'bold/009';
'bold/018';
'bold/019';
'bold/020';
'bold/021';
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for d_idx=1:length(dirs)

	mc=load(sprintf('%s/%s/%s',bold_dir,dirs{d_idx},fmc_stem));

	mc_regressor{d_idx}=mc(:,2:7);

end;
mc_dirs=dirs;

save mc_regressor.mat mc_regressor mc_dirs

fprintf('DONE!\n');

return;
