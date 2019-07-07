close all; clear all;

%------------------------------------------------------------------------------------
%----------------------------------GLM setup-----------------------------------------
% er-fmri setup

file_stc={
    '063019_2_fsaverage_32ch_mb_run_1';
    '063019_2_fsaverage_32ch_mb_run_2';
    };


TR=2.0; %second

exclude_time=[1:2];
%exclude_time=[];

confound_polynomial_order=2;

output_stem='tsnr_32ch';
%------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% construct contrast matrix using legendres to model the effects
%

fprintf('loading data...\n');
stc=[];
exclude_time_all=[];
for f_idx=1:length(file_stc)
    fn=sprintf('%s-lh.stc',file_stc{f_idx});
    fprintf('<%s>...\r',fn);
    [stc_lh,v_lh,cc,dd]=inverse_read_stc(fn);
    fn=sprintf('%s-rh.stc',file_stc{f_idx});
    fprintf('<%s>...\r',fn);
    [stc_rh,v_rh,cc,dd]=inverse_read_stc(fn);
    
    stc=cat(2,stc,cat(1,stc_lh,stc_rh));
    if(~isempty(exclude_time))
        exclude_time_all=cat(1,exclude_time_all,exclude_time(:)+(size(stc_lh,2)).*(f_idx-1));
    end;
    
    timepoints(f_idx)=size(stc_lh,2);
    
end;
stc(:,exclude_time_all)=[];
fprintf('\n');
stc=stc';

idx=find(isnan(stc(:)));
stc(idx)=randn(size(idx)).*eps;

idx=find(abs(stc(:))<eps);
stc(idx)=randn(size(idx)).*eps;

cumsum_timepoints=cumsum(timepoints);
cumsum_timepoints=cat(1,0,cumsum_timepoints(:));

fprintf('adding confounds into a contrast matrix...\n');
n_confound=1;
confound=[];
beta_dc=[];
confound_period=[60, 132, 180]; %second
confound_period=[];


contrast_count=0;

for run_idx=1:length(file_stc)
    for j=0:confound_polynomial_order
        contrast_count=contrast_count+1;
        
        if(j==0) beta_dc=cat(1,beta_dc,contrast_count); end;
        
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound)=([0:1/((timepoints(run_idx))-1):1].^(j))';
        n_confound=n_confound+1;
    end;
    timeVec=[0:TR:TR*(timepoints-1)];
    for j=1:length(confound_period)
        contrast_count=contrast_count+1;
        contrast_count=contrast_count+1;
        
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound)=cos(2.*pi./confound_period(j).*timeVec)';
        n_confound=n_confound+1;
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound)=sin(2.*pi./confound_period(j).*timeVec)';
        n_confound=n_confound+1;
    end;
end;

%remove time points
confound(exclude_time_all,:)=[];

% %global mean
% confound(:,end+1)=mean(stc,2);
% n_confound=n_confound+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove confound
beta=inv(confound'*confound)*confound'*stc;
res=stc-confound*inv(confound'*confound)*confound'*stc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg_idx=[1:confound_polynomial_order+1:size(confound,2)];
avg=mean(beta(avg_idx,:));

tsnr=avg./std(res,0,1);

fprintf('\tarchiving results...\n');
fn=sprintf('%s_tsnr-lh.stc',output_stem);
inverse_write_stc(repmat(tsnr(:,1:length(v_lh))',[1 5]),v_lh,0,TR.*1e3,fn);

fn=sprintf('%s_tsnr-rh.stc',output_stem);
inverse_write_stc(repmat(tsnr(:,length(v_lh)+1:end)',[1 5]),v_rh,0,TR.*1e3,fn);
return;
    

