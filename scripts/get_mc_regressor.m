close all; clear all;

fmc_stem='fmcpr.mcdat';

bold_dir='/space/maki7/1/fhlin/seeg/s026/fmri_data/unpack';

dirs={
'bold/030';
'bold/032';
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for d_idx=1:length(dirs)

	mc=load(sprintf('%s/%s/%s',bold_dir,dirs{d_idx},fmc_stem));

	mc_regressor{d_idx}=mc(:,2:7);

end;
save mc_regressor.mat mc_regressor

fprintf('DONE!\n');

return;
