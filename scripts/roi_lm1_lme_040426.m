close all; clear all;

% --- 1. Load the Data ---
filename = 'roi_data_lm1.xlsx';

% Preserve original column names so '10-Hz-3-ppb' doesn't get corrupted by MATLAB
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
T = readtable(filename, opts);

% Extract subjects and data
subjects_raw = T{:, 1};       % First column is subject IDs
data_matrix = T{:, 2:6};      % Columns 2 through 6 are the ROI values

% The 6 condition names exactly as they appear in the file header
cond_str = T.Properties.VariableNames(2:6)'; 

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


C_str = {
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