close all; clear all;

%------------------------------------------------------------------------------------
%----------------------------------GLM setup-----------------------------------------
% er-fmri setup

file_stc={
    '../fmri_data/unpack/bold/030/s026_2_fsaverage_sfmcprstc';
    '../fmri_data/unpack/bold/032/s026_2_fsaverage_sfmcprstc';
    };


hdr_length=30.0; %second
hdr_pre=6.0; 	%second (pre-stim)
TR=2; %second
hdr_timeVec=[0:TR:hdr_length-TR]-hdr_pre;
n_hdr=round(hdr_length/TR);


flag_hdr_canonical=1;
flag_hdr_fir=0;

erfmri_para={
    'fmri_soa_01.para';
    'fmri_soa_02.para';
    };

file_mc_regressor='mc_regressor.mat'; %motion correction coefficients for regession

file_ventrical_wm={
    'regressor_wm_ventrical_030.mat';
    'regressor_wm_ventrical_032.mat';
    };


%describe the variables to be modeled in GLM; one cell element represents
%one HDR type to be estimated
glm_var={
    [1];
    [2];
    [3];
    };

n_run=length(erfmri_para);

%auto-regressive modeling of fMRI time series residual
flag_ar=0;
flag_ar_global=1;
max_ar_order=1;

%exclude_time=[1:2];
exclude_time=[];

confound_polynomial_order=0;

hypothesis{1}.rv=[1];	
hypothesis{1}.cvec=[1];

hypothesis{2}.rv=[2];	
hypothesis{2}.cvec=[2];

hypothesis{3}.rv=[3];	
hypothesis{3}.cvec=[3];


output_stem='fmri_surf_soa_glm';
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
        exclude_time_all=cat(1,exclude_time_all,exclude_time(:)+(size(stc0,2)).*(f_idx-1));
    end;
    
    timepoints(f_idx)=size(stc_lh,2);
    
end;
stc(:,exclude_time_all)=[];
fprintf('\n');
stc=stc';

idx=find(isnan(stc(:)));
stc(idx)=randn(size(idx)).*eps;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% construct contrast matrix using legendres to model the effects
%

%loading behavior data
fprintf('preparation of contrast matrix...\n');

para_time=[];
para_type=[];
cumsum_timepoints=cumsum(timepoints);
cumsum_timepoints=cat(1,0,cumsum_timepoints(:));

for para_count=1:length(erfmri_para)
    fprintf('loading paradigm SOA [%s]...\n',erfmri_para{para_count});
    [p_t,p_p]=textread(erfmri_para{para_count},'%f %d');
    para_time=[para_time; p_t+TR*cumsum_timepoints(para_count)];
    para_type=[para_type; p_p];
end;


fprintf('analyzing paradigms and creating contrast matrix...\n');

para_all=unique(para_type);
scm=[];
for glm_var_idx=1:length(glm_var)
    glm_var_now=glm_var{glm_var_idx};
    
    idx=[];
    for g_idx=1:length(glm_var_now)
        idx=union(idx,find(para_type==glm_var_now(g_idx)));
    end;
    if(~isempty(idx))
        fprintf('paradigm type [%s] found! [%d] occurances...\n',mat2str(glm_var_now),length(idx));
        if(isempty(scm))
            scm=zeros(sum(timepoints)*20,1); %20-fold temporal over-sampling
        end;
        
        soa_idx=round(para_time(idx)/TR*20);
        soa_idx(find(soa_idx<=0))=[];
        soa_idx(find(soa_idx>size(scm,1)))=[];
        scm(soa_idx,glm_var_idx)=1;
    end;
end;


%FIR HDR
if(flag_hdr_canonical)
    hhdr_timeVec=[0:TR/20:hdr_length-TR]; %20-fold over-sampling
    HDR=fmri_hdr(hhdr_timeVec,'hdr_order',1);
    HDR(find(hhdr_timeVec<0),:)=0;
    tmp=sum(HDR.^2,1);
    HDR=HDR./repmat(sqrt(tmp),[size(HDR,1),1]);
    
    hhdr_timeVec0=[0:TR:hdr_length-TR];
    HDR0=fmri_hdr(hhdr_timeVec0,'hdr_order',1);
    HDR0(find(hhdr_timeVec0<0),:)=0;
    tmp=sum(HDR0.^2,1);
    HDR0=HDR0./repmat(sqrt(tmp),[size(HDR0,1),1]);
elseif(flag_hdr_fir)
    HDR=kron(eye(n_hdr),ones(20,1));
    HDR0=eye(n_hdr);
end;


contrast_hdr=[];
for c_idx=1:size(scm,2)
    for hdr_idx=1:size(HDR,2)
        onset=scm(:,c_idx);
        if(flag_hdr_fir)
            onset=circshift(onset,-round(hdr_pre/TR*20)+1);
        end;
        tmp=conv(onset,HDR(:,hdr_idx));
        tmp=tmp(1:size(scm,1));
        contrast_hdr(:,end+1)=tmp(:);
    end;
end;

if(flag_hdr_canonical)
    contrast_hdr=contrast_hdr(1:20:end,:);
elseif(flag_hdr_fir)
    contrast_hdr=contrast_hdr(1:20:end,:);
end;
contrast_hdr(exclude_time_all,:)=[];

contrast=contrast_hdr;
contrast_of_interest=size(scm,2);

contrast_count=size(contrast_hdr,2);
fprintf('adding confounds into contrast matrix...\n');
n_confound=1;
confound=[];
beta_dc=[];
confound_period=[60, 120, 180]; %second
confound_period=[];


%motion regressors
if(~isempty(file_mc_regressor))
    load(file_mc_regressor);
    for run_idx=1:n_run
        mcr=mc_regressor{run_idx};
        mcr(end,:)=[];
        %mcr([exclude_time(:); exclude_time(end)+1],:)=[]; %one more dummy because stc files have one-less time points
        %mcr([exclude_time(:)],:)=[]; %one more dummy because stc files have one-less time points
        
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound:n_confound+5)=fmri_scale(mcr,1,0);
        n_confound=n_confound+6;
    end;
end;


%white matter/ventricle regressors
if(~isempty(file_ventrical_wm))
    for run_idx=1:n_run
        load(file_ventrical_wm{run_idx});
        mcr=regressor_ventricle(:);    
        mcr(end,:)=[];
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound:n_confound)=fmri_scale(mcr,1,0);
        
        mcr=regressor_wm(:);    
        mcr(end,:)=[];
        confound(1+cumsum_timepoints(run_idx):cumsum_timepoints(run_idx+1),n_confound+1:n_confound+1)=fmri_scale(mcr,1,0);

        n_confound=n_confound+2;
    end;
end;

for run_idx=1:n_run
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

%global mean
%confound(:,end+1)=mean(stc,2);
%n_confound=n_confound+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%remove confound
%stc=stc-confound*inv(confound'*confound)*confound'*stc;
%error_sig02=sum(stc.^2,1)./(size(stc,1)-size(confound,2));
contrast=cat(2,contrast,confound);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('GLM estimation....');
beta=inv(contrast'*contrast)*contrast'*stc;
error=(stc-contrast*beta);
error_sig2=sum(error.^2,1)./(size(stc,1)-size(contrast,2));
df=size(stc,1)-rank(contrast);


if(flag_ar)
    %estimate noise auto-correlation
    %h = waitbar(0,'auto-correlation on residual and re-estimating...');
    %warning off;
    if(flag_ar_global)
        ii=randsample(size(error,1),100);
        ee=error(ii,:)';
        [w,a,c]=arfit(ee(:),1,max_ar_order);
        rho=[1 a];
        
        [Ainvt postdef]=chol(toeplitz(rho));
        p1=size(Ainvt,1);
        A=inv(Ainvt');
        Vmhalf=toeplitz([A(p1, p1:-1:1) zeros(1,size(error,1)-p1)],zeros(1,size(error,1)));
        Vmhalf(1:p1,1:p1)=A;
        
        %generate whitened measuremends and design matrix; OLS estimation
        w_contrast=Vmhalf*contrast;
        w_y=Vmhalf*stc;
        w_beta=inv(w_contrast'*w_contrast)*w_contrast'*w_y;
        w_dof=size(error,1)-rank(w_contrast);
        w_res=w_y-w_contrast*w_beta;
        w_error_sig2=sum(w_res.^2,1)./w_dof;
    else        
        for v_idx=1:size(error,2)
            %waitbar(v_idx/size(error,2),h);
            if(mod(v_idx,1)==0) fprintf('auto-correlation estimation :: [%05d|%05d::%1.1f%%]...\r',v_idx,size(error,2),v_idx/size(error,2)*100); end;
            %[acf(:,v_idx),lags] = autocorr(error(:,v_idx) , 10);
            [w,a,c]=arfit(error(:,v_idx),1,max_ar_order);
            rho=[1 a];
            
            [Ainvt postdef]=chol(toeplitz(rho));
            p1=size(Ainvt,1);
            A=inv(Ainvt');
            Vmhalf=toeplitz([A(p1, p1:-1:1) zeros(1,size(error,1)-p1)],zeros(1,size(error,1)));
            Vmhalf(1:p1,1:p1)=A;
            
            %generate whitened measuremends and design matrix; OLS estimation
            w_contrast=Vmhalf*contrast;
            w_y=Vmhalf*stc(:,v_idx);
            w_beta(:,v_idx)=inv(w_contrast'*w_contrast)*w_contrast'*w_y;
            w_dof(v_idx)=size(error,1)-rank(w_contrast);
            w_res=w_y-w_contrast*w_beta(:,v_idx);
            w_error_sig2(v_idx)=(w_res'*w_res)/w_dof(v_idx);
        end;
    end;
    %close(h);
else

    w_beta=beta;
    w_error_sig2=error_sig2;
    w_contrast=contrast;

end;


% if(size(contrast,2)>1)%F-test only when rank(contrast)>1; t-test suffices otherwise.
%     f_stat=(error_sig02-error_sig2)./error_sig2.*(size(contrast,1)-size(contrast,2))./(size(contrast,2)-1);
%     fn=sprintf('%s_fstat-lh.stc',output_stem);
%     inverse_write_stc(repmat(f_stat(1:length(v_lh))',[1 5]),v_lh,0,TR.*1e3,fn);
%     fn=sprintf('%s_fstat-rh.stc',output_stem);
%     inverse_write_stc(repmat(f_stat(length(v_lh)+1:end)',[1 5]),v_rh,0,TR.*1e3,fn);
%     
%     figure;
%     h=etc_render_fsbrain_stc({sprintf('%s_fstat',output_stem)},[25 40]);
% end;
fprintf('statistical hypothesis test...\n');
for h=1:length(hypothesis)
    
    fprintf('hypothesis testing %s with weights %s...\n',mat2str(glm_var{hypothesis{h}.rv}),mat2str(hypothesis{h}.cvec));
    
    %build up contrast vector
    if(~isempty(hypothesis{h}.rv))
                rv_idx=hypothesis{h}.rv;
                fprintf('\t %s\n', mat2str(glm_var{rv_idx}));
    end;
    
    cc_hdr=zeros(size(scm,2),1);
    %cc_hdr(hypothesis{h}.rv)=hypothesis{h}.cvec;
    ii=[];
    for rv_idx=1:length(hypothesis{h}.rv)
        ii(rv_idx) = find([glm_var{:}] == hypothesis{h}.rv(rv_idx));
    end;
    cc_hdr(ii)=hypothesis{h}.cvec;
    
    
    cc=zeros(size(w_beta,1),1);
    cc_tmp=kron(cc_hdr(:),ones(size(HDR0,2),1));
    cc(1:length(cc_tmp(:)))=cc_tmp(:);
    t_stat=cc'*w_beta./sqrt(w_error_sig2.*(cc'*inv(w_contrast'*w_contrast)*cc));
    
    
    fprintf('\tarchiving results...\n');
    %t_stat=reshape(t_stat,[nii_sz(1) nii_sz(2) nii_sz(3)]);
    %etc_save_nii(t_stat,nii,sprintf('%s_tstat.nii',output_stem));

    %beta=reshape(beta,[nii_sz(1) nii_sz(2) nii_sz(3) size(beta,1)]);
    %etc_save_nii(beta,nii,sprintf('%s_beta.nii',output_stem));
    
    fn=sprintf('%s_h%02d_tstat-lh.stc',output_stem,h);
    inverse_write_stc(repmat(t_stat(1:length(v_lh))',[1 5]),v_lh,0,TR.*1e3,fn);
    fn=sprintf('%s_h%02d_tstat-rh.stc',output_stem,h);
    inverse_write_stc(repmat(t_stat(length(v_lh)+1:end)',[1 5]),v_rh,0,TR.*1e3,fn);
    
    fn=sprintf('%s_h%02d_beta-lh.stc',output_stem,h);
    inverse_write_stc(beta(:,1:length(v_lh))',v_lh,0,TR.*1e3,fn);
    fn=sprintf('%s_h%02d_beta-rh.stc',output_stem,h);
    inverse_write_stc(beta(:,length(v_lh)+1:end)',v_rh,0,TR.*1e3,fn);
    
    if(flag_hdr_fir)
        dspm=zeros(length(hdr_timeVec),size(stc,2));
        mne=zeros(length(hdr_timeVec),size(stc,2));
        
        fprintf('\tdynamic SPM for FIR HDR analysis...\n');
        for t_idx=1:length(hdr_timeVec)
            fprintf('*');
            
            cc=zeros(size(w_beta,1),1);
            tmp=zeros(size(HDR0,2),1);
            tmp(t_idx)=1;
            cc_tmp=kron(cc_hdr(:),tmp(:));
            cc(1:length(cc_tmp(:)))=cc_tmp(:);
            
            mne(t_idx,:)=cc'*w_beta;
        end;
        baseline_idx=find(hdr_timeVec<0);
        dspm=mne./repmat(std(mne(baseline_idx,:),[],1),[size(mne,1),1]);
        fprintf('\n');
        
        fn=sprintf('%s_h%02d_mne-lh.stc',output_stem,h);
        inverse_write_stc(mne(:,1:length(v_lh))',v_lh,hdr_timeVec(1).*1e3,TR.*1e3,fn);
        fn=sprintf('%s_h%02d_mne-rh.stc',output_stem,h);
        inverse_write_stc(mne(:,length(v_lh)+1:end)',v_rh,hdr_timeVec(1).*1e3,TR.*1e3,fn);
        
        fn=sprintf('%s_h%02d_dspm-lh.stc',output_stem,h);
        inverse_write_stc(dspm(:,1:length(v_lh))',v_lh,hdr_timeVec(1).*1e3,TR.*1e3,fn);
        fn=sprintf('%s_h%02d_dspm-rh.stc',output_stem,h);
        inverse_write_stc(dspm(:,length(v_lh)+1:end)',v_rh,hdr_timeVec(1).*1e3,TR.*1e3,fn);
        
    end;
end;

fprintf('DONE\n');

return;


