close all; clear all;

fmc_stem='fmcpr.mcdat';

bold_dir='/space_lin2/fhlin/eegmri_memory/s012/resting_data/unpack';
bold_dir='/Users/fhlin/workspace/eegmri_memory/s012/resting_data/unpack';

dirs={
'bold/005';
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for d_idx=1:length(dirs)

	mc=load(sprintf('%s/%s/%s',bold_dir,dirs{d_idx},fmc_stem));

	mc_regressor{d_idx}=mc(:,2:7);

        mc_regressor_this_run=mc(:,2:7);

        idx=findstr(dirs{d_idx},'/');
        runstr=dirs{d_idx}(idx+1:end);
        save(sprintf('mc_regressor_%s.mat',runstr),'mc_regressor_this_run');
end;
mc_dirs=dirs;

save mc_regressor.mat mc_regressor mc_dirs

fprintf('DONE!\n');

return;
