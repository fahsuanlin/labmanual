close all; clear all;

% --- 1. Load the Data ---
filename = 'meg_n8.xlsx';

% Preserve original column names so '10-Hz-3-ppb' doesn't get corrupted by MATLAB
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
T = readtable(filename, opts);

% Extract subjects and data
subjects_raw = T{:, 1};       % First column is subject IDs
data_matrix = T{:, 2:7};      % Columns 2 through 7 are the MEG values

% The 6 condition names exactly as they appear in the file header
cond_str = T.Properties.VariableNames(2:7)'; 

n_sub = length(subjects_raw);
n_cond = length(cond_str);

% --- 2. Reshape into Long Format for LME ---
% Flatten the data matrix column by column
y = data_matrix(:); 

% Replicate subject labels (Sub1..8, Sub1..8...)
Subjects = categorical(repmat(subjects_raw, n_cond, 1));

% Replicate condition labels (Baselinex8, 0.5-Hzx8...)
tmp_cond = repelem(cond_str, n_sub);
% Explicitly set categorical order so Baseline is locked as the Intercept
Conditions = categorical(tmp_cond, cond_str); 

% Build the Master Table
tbl = table(y, Conditions, Subjects, 'VariableNames', {'Con', 'Cond', 'Subj'});

% --- 3. Fit the LME Model ---
fprintf('Fitting Linear Mixed Effects Model on %d subjects...\n', n_sub);
warning('off', 'stats:LinearMixedModel:FitWarning');
lme = fitlme(tbl, 'Con ~ Cond + (1|Subj)', 'FitMethod', 'REML');
warning('on', 'stats:LinearMixedModel:FitWarning');

% Display the raw model
disp(lme);

% --- 4. Define Contrasts Matrix (6 Columns Now) ---
% Order: [Intercept(Baseline), Diff_0.5Hz, Diff_10Hz3, Diff_10Hz5, Diff_30Hz3, Diff_30Hz5]
C_matrix = [
    1  0  0  0  0  0;  % C1: Baseline Mean
    1  1  0  0  0  0;  % C2: 0.5-Hz Mean
    1  0  1  0  0  0;  % C3: 10-Hz-3-ppb Mean
    1  0  0  1  0  0;  % C4: 10-Hz-5-ppb Mean
    1  0  0  0  1  0;  % C5: 30-Hz-3-ppb Mean
    1  0  0  0  0  1;  % C6: 30-Hz-5-ppb Mean
    0  1  0  0  0  0;  % C7: 0.5-Hz vs Baseline
    0 -1  1  0  0  0;  % C8: 10-Hz-3-ppb vs 0.5-Hz  (Beta 3 - Beta 2)
    0 -1  0  1  0  0;  % C9: 10-Hz-5-ppb vs 0.5-Hz  (Beta 4 - Beta 2)
    0 -1  0  0  1  0;  % C10: 30-Hz-3-ppb vs 0.5-Hz (Beta 5 - Beta 2)
    0 -1  0  0  0  1;  % C11: 30-Hz-5-ppb vs 0.5-Hz (Beta 6 - Beta 2)
];

C_str = {
    'Baseline Mean';
    '0.5-Hz Mean';
    '10Hz-3ppb Mean';
    '10Hz-5ppb Mean';
    '30Hz-3ppb Mean';
    '30Hz-5ppb Mean';
    '0.5-Hz vs Baseline';
    '10Hz-3ppb vs 0.5-Hz';
    '10Hz-5ppb vs 0.5-Hz';
    '30Hz-3ppb vs 0.5-Hz';
    '30Hz-5ppb vs 0.5-Hz';
};

n_contrasts = size(C_matrix, 1);

% --- 5. Evaluate and Print Contrasts ---
fprintf('\n================== CONTRAST RESULTS ==================\n');
betas = fixedEffects(lme);

for c_idx = 1:n_contrasts
    C = C_matrix(c_idx, :);
    
    % Calculate the exact effect size for this contrast
    est = sum(C(:) .* betas(:)); 
    
    % Run Satterthwaite test
    [p_val, F_stat, ~] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
    
    % Recover the true t-statistic
    t_stat = sign(est) * sqrt(F_stat);
    
    % Print neatly to the console
    fprintf('Contrast %2d (%-20s): Effect = %8.3f | t = %7.3f | p = %.5f', ...
        c_idx, C_str{c_idx}, est, t_stat, p_val);
        
    % Add a visual star if it's significant (p < 0.05)
    if p_val < 0.05
        fprintf('  ***\n');
    else
        fprintf('\n');
    end
end
fprintf('======================================================\n');