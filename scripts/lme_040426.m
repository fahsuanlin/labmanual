close all; clear all;

subject={
    'sub-101';
    'sub-102';
    'sub-103';
    'sub-104';
    'sub-106';
    'sub-107';
    'sub-109';
    'sub-111';
    'sub-112';
    };

cond_nii={
    'con_0002.nii';
    'con_0003.nii';
    'con_0004.nii';
    'con_0005.nii';
    'con_0006.nii';
    };

cond_str={
    '0.5-Hz';
    '10-Hz-3-ppb';
    '10-Hz-5-ppb';
    '30-Hz-3-ppb';
    '30-Hz-5-ppb';
    };


% --- 1. Define your Contrasts Matrix (UPDATED FOR INTERCEPT MODEL) ---
% Model order: [Intercept(0.5Hz), Diff_10Hz3, Diff_10Hz5, Diff_30Hz3, Diff_30Hz5]
C_matrix = [
    1  0    0    0    0;    % Contrast 1: 0.5Hz Mean (Intercept)
    1  1    0    0    0;    % Contrast 2: 10Hz-3ppb Mean (Int + Diff2)
    1  0    1    0    0;    % Contrast 3: 10Hz-5ppb Mean (Int + Diff3)
    1  0    0    1    0;    % Contrast 4: 30Hz-3ppb Mean (Int + Diff4)
    1  0    0    0    1;    % Contrast 5: 30Hz-5ppb Mean (Int + Diff5)
    0  0.5  0.5  0    0;    % Contrast 6: 10Hz vs 0.5Hz 
    0  0    0    0.5  0.5;  % Contrast 7: 30Hz vs 0.5Hz
    0  0.25 0.25 0.25 0.25; % Contrast 8: 10Hz/30Hz vs 0.5Hz
    0 -1   -1    1    1;    % Contrast 9: 30Hz vs 10Hz
];


C_str={
    '0.5hz';
    '10hz-3ppb';
    '10hz-5ppb';
    '30hz-3ppb';
    '30hz-5ppb';
    '10hz-0.5hz';
    '30hz-0.5hz';
    '10hz30hz-0.5hz';
    '30hz-10hz';
    };


stc=[];
for subj_idx=1:length(subject)
    for cond_idx=1:length(cond_nii)
        data=MRIread(sprintf('%s/%s',subject{subj_idx},cond_nii{cond_idx}));

        stc=cat(1,stc,data.vol(:)');

    end;
end;
n_voxels=size(stc,2);
valid_voxels=find(~isnan(sum(stc,1)));

tmp=repmat(subject(:),1,5)';
Subjects = categorical(tmp(:));

tmp=repmat(cond_str(:),1,9);
Conditions = categorical(tmp(:));


n_contrasts = size(C_matrix, 1);

% --- 2. Preallocate Arrays ---
tstat_lme = zeros(n_voxels, n_contrasts);
effect_lme = zeros(n_voxels, n_contrasts);
pval_lme = ones(n_voxels, n_contrasts); 
qval_lme = ones(n_voxels, n_contrasts); 

% Temporary sliced arrays strictly for the parfor loop
num_valid = length(valid_voxels);
par_tstat_all = zeros(num_valid, n_contrasts);
par_effect_all = zeros(num_valid, n_contrasts);
par_pval_all = ones(num_valid, n_contrasts); 

fprintf('Running parallel LME on %d valid voxels for %d contrasts...\n', num_valid, n_contrasts);

% --- 3. Parallel LME Loop ---
parfor i = 1:num_valid
    
    v_idx = valid_voxels(i);
    y = stc(:, v_idx);
    
    % Skip flatline voxels outside the brain mask
    if var(y) == 0
        continue; 
    end
    
    % Build table
    tbl = table(y, Conditions, Subjects, 'VariableNames', {'Con', 'Cond', 'Subj'});
    
    % Fit the DEFAULT Intercept model (Guaranteed 5 coefficients)
    warning('off', 'stats:LinearMixedModel:FitWarning');
    lme = fitlme(tbl, 'Con ~ Cond + (1|Subj)', 'FitMethod', 'REML');
    warning('on', 'stats:LinearMixedModel:FitWarning');
    
    betas = fixedEffects(lme); 
    
    par_tstat = zeros(1, n_contrasts);
    par_effect = zeros(1, n_contrasts);
    par_pval = ones(1, n_contrasts);
    
    % --- Evaluate ALL Contrasts ---
    for c_idx = 1:n_contrasts
        C = C_matrix(c_idx, :);
        
        % Calculate Effect Size
        est = sum(C(:) .* betas(:)); 
        
        % Satterthwaite test with a safety fallback
        try
            [p_val, F_stat, ~] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
        catch
            % If variance is too low for Satterthwaite, fallback to residual DF
            [p_val, F_stat, ~] = coefTest(lme, C, 0, 'DFMethod', 'residual');
        end
        
        % Recover the true t-statistic
        t_stat = sign(est) * sqrt(F_stat);
        
        par_tstat(c_idx) = t_stat;
        par_effect(c_idx) = est;
        par_pval(c_idx) = p_val;
    end
    
    par_tstat_all(i, :) = par_tstat;
    par_effect_all(i, :) = par_effect;
    par_pval_all(i, :) = par_pval;
end

% --- 4. Reconstruct Full Brain Maps ---
tstat_lme(valid_voxels, :) = par_tstat_all;
effect_lme(valid_voxels, :) = par_effect_all;
pval_lme(valid_voxels, :) = par_pval_all;


for c_idx=1:n_contrasts
    [pp,qq]=FDR(par_pval_all(:,c_idx),0.05);
    if(~isempty(pp))
        % 1. Find indices of all voxels that survived the threshold
        sig_idx = find(pval_lme(:, c_idx) <= pp);
        num_sig = length(sig_idx);

        % 2. Extract their t-statistics to find the empirical boundary
        surviving_tstats = tstat_lme(sig_idx, c_idx);
        t_threshold = min(abs(surviving_tstats)); % The lowest |t| that survived

        fprintf('Contrast [%d] (%s): %d voxels survive (p < %g | |t| > %.3f)\n', ...
            c_idx, C_str{c_idx}, num_sig, pp, t_threshold);

        % 3. Create a strict FDR-Masked map for visualization
        masked_tstat = zeros(n_voxels, 1);
        masked_tstat(sig_idx) = tstat_lme(sig_idx, c_idx); % Copy ONLY significant t-stats

        data.vol = reshape(masked_tstat, data.volsize);
        data.nframes = 1;
        MRIwrite(data, sprintf('%s_tstat_FDR05_masked.nii', C_str{c_idx}));

    else
        fprintf('Contrast [%d] (%s): 0 voxels survive FDR q<0.05\n', c_idx, C_str{c_idx});

        % Optionally output a completely blank map for consistency
        data.vol = zeros(data.volsize);
        data.nframes = 1;
        MRIwrite(data, sprintf('%s_tstat_FDR05_masked.nii', C_str{c_idx}));
    end
end;


% --- 6. Export Results ---
fprintf('Exporting NIfTI files...\n');
for c_idx = 1:n_contrasts
    % Export T-statistics
    data.vol = reshape(tstat_lme(:,c_idx), data.volsize);
    data.nframes = 1;
    MRIwrite(data, sprintf('%s_tstat.nii', C_str{c_idx}));
    
    % Export P-values
    data.vol = reshape(pval_lme(:,c_idx), data.volsize);
    MRIwrite(data, sprintf('%s_pval.nii', C_str{c_idx}));

    % Export Q-values
    data.vol = reshape(qval_lme(:,c_idx), data.volsize);
    MRIwrite(data, sprintf('%s_qval.nii', C_str{c_idx}));
end

% Save everything to the .mat file
save lme_040426 tstat_lme pval_lme qval_lme effect_lme C_matrix C_str;