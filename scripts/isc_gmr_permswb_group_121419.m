close all; clear all;

subj_path={
% %   '180810_social';
% %   '181008_social';
    '181120_social';
    '190109_social';
%    '190116_social';
    '190215_social';
    '190315_social';
    '190426_social';
    '190524_social';
    '190531_social';
    '190621_social';
    '190726_social';
%    '190823_social';
    '190830_social';
    '190906_social';
    '190920_social';
    '190927_social';
%    '191004_social';
    '191018_social';
    '191019_social';
    '191101_social';
    '191115_social';
    };

% BL: 1
% BH: 2
% CL: 1
% CH: 2
% GL: 1
% GH: 2

movie_sequence={
    [1 4 5 3 6 2]; %181120_social
    [3 6 1 5 2 4]; %190109_social
%    [6 1 4 2 3 5]; %190116_social
    [2 3 6 4 5 1]; %190215_social
    [4 2 5 1 6 3]; %190315_social
    [5 4 1 2 3 6]; %190426_social
    [3 1 6 5 4 2]; %190524_social
    [6 3 2 4 1 5]; %190531_social
    [6 1 3 2 5 4]; %190621_social
    [5 3 2 1 6 4]; %190726_social
%    [4 1 5 6 3 2]; %190823_social
    [1 5 3 2 4 6]; %190830_social
    [3 2 6 1 4 5]; %190906_social
    [6 4 1 5 2 3]; %190920_social
    [1 4 6 3 5 2]; %190927_social
%    [3 6 2 5 1 4]; %191004_social
    [5 2 4 1 3 6]; %191018_social
    [2 5 3 6 4 1]; %191019_social
    [4 1 5 2 6 3]; %191101_social
    [6 3 1 4 2 5]; %191115_social
    };



TR=2;


hemi={
    'lh';
    'rh';
    };

sequence=[1 2 3 4 5 6]; %three different movie sequences shown to subjects; ISC calculated separately
 
hypothesis={
    [1 1 1 1 1 1]; %grand average
    [1 0 1 0 1 0]; %low social
    [0 1 0 1 0 1]; %high social
    [-1 1 -1 1 -1 1]; %high - low social
    };


root_path='/Users/fhlin_admin/workspace/social_eegmri/';

output_stem='isc_permswb_group_121419_z';

confound_polynomial_order=2;
confound_sinusoidal_order=3;

for sequence_idx=1:length(sequence)
    %subject=subj_path(find(movie_sequence==sequence(sequence_idx)));
    
    
    %oo=sprintf('%s_sequence%03d_%s',output_stem,sequence(sequence_idx),run_stem);
    oo=sprintf('%s_sequence%02d',output_stem,sequence(sequence_idx));
    subject=subj_path;
    
    for hemi_idx=1:length(hemi)
        stc=[];
        for subj_idx=1:length(subject)
            
            run_idx=find(movie_sequence{subj_idx}==sequence(sequence_idx));
            
            [stc(:,:,subj_idx),a,b,c]=inverse_read_stc(sprintf('%s/%s/fmri_analysis/%s_2_fsaverage_mb_run_%d_acc-%s.stc',root_path,subject{subj_idx},subject{subj_idx},run_idx,hemi{hemi_idx}));
            fprintf('<%s>::<%s>\t%05d time ponits\n',subject{subj_idx},hemi{hemi_idx},size(stc,2));
            %			length(find(isnan(stc(:))))
        end;
        ii=find(isnan(stc(:)));
        stc(ii)=randn(size(ii)).*eps;
        
        fprintf('\n');
        
        %remove the first 5 time points
        %stc=stc(:,6:end,:);\
        stc=stc(:,11:end,:);
        %stc(:,film_set_fixation{sequence_idx},:)=[]; %remove inter-film blank interval and first 15 s.
        
        %remove global mean
        for s_idx=1:size(stc,3)
            tmp=squeeze(stc(:,:,s_idx));
            
            ga=mean(tmp,1);
            
            tmp=tmp-tmp*ga'*inv(ga*ga')*ga;
            
            stc(:,:,s_idx)=tmp;
        end;
        
        timeVec=[1:size(stc,2)]';
        D_poly=[];
        D_sinu=[];
        D_poly=ones(length(timeVec),1);
        for i=1:confound_polynomial_order
            tmp=timeVec.^(i);
            D_poly(:,i+1)=fmri_scale(tmp(:),1,0);
        end;
        for i=1:confound_sinusoidal_order
            D_sinu(:,i*2-1)=sin(timeVec.*i./timeVec(end).*pi);
            D_sinu(:,i*2)=cos(timeVec.*i./timeVec(end).*pi);
        end;
        D=cat(2,D_poly,D_sinu);
        
        if(~isempty(D))
            D_prep=D*inv(D'*D)*D';
        else
            D_prep=[];
        end;
        
        for v_idx=1:size(stc,1)
            if(mod(v_idx,100)==0)
                fprintf('[%1.1f%%]...\r',v_idx/size(stc,1)*100);
            end;
            if(~isempty(D_prep))
                data=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:));
            else
                data=squeeze(stc(v_idx,:,:));
            end;
            ii=find(abs(data)<eps);
            data(ii)=randn(size(ii)).*eps;
            
            %save all correlation coef.
            if((sequence_idx==1)&&(v_idx==1))
                corr_buffer{hemi_idx}=zeros(length(sequence),size(stc,1),size(stc,3),size(stc,3));
            end;
            corr_buffer{hemi_idx}(sequence_idx,v_idx,:,:)=corrcoef(data);
            
            tmp=tril(corrcoef(data),-1);
            rr=tmp(find(tmp));
            z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
            zm=mean(z);
            zmd=median(z);

            if((sequence_idx==1)&&(v_idx==1))
                cc{hemi_idx}=zeros(length(sequence),size(stc,1),length(z(:)));
            end;
            cc{hemi_idx}(sequence_idx,v_idx,:)=z(:)';
            if(find(isnan(z))) keyboard; end;
            
            ccm(v_idx)=zm;
            ccmd(v_idx)=zmd;
            
        end;
        fprintf('\n');
        
        %if(~isempty(find(isnan(cc{sequence_idx}(:))))) keyboard; end;
        
        inverse_write_stc(squeeze(cc{hemi_idx}(sequence_idx,:,:)),a,b,c,sprintf('%s_gmr-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(ccm(:),[1 5]),a,b,c,sprintf('%s_gmr_mean-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(ccmd(:),[1 5]),a,b,c,sprintf('%s_gmr_median-%s.stc',oo,hemi{hemi_idx}));
    end;
end;

for hypothesis_idx=1:length(hypothesis)
    for hemi_idx=1:length(hemi)        %permutation
        for v_idx=1:size(stc,1)
            tmp=hypothesis{hypothesis_idx}(:)'*squeeze(cc{hemi_idx}(:,v_idx,:));
            effect_mean{hemi_idx}(hypothesis_idx,v_idx)=mean(tmp);
            effect_median{hemi_idx}(hypothesis_idx,v_idx)=median(tmp);
        end;
    end;
end;


for hypothesis_idx=1:length(hypothesis)
    
    oo=sprintf('%s_hypothesis%02d',output_stem,hypothesis_idx);

    for hemi_idx=1:length(hemi)        %permutation
        zm_perm_p=zeros(size(ccm));
        zmd_perm_p=zeros(size(ccmd));
        n_perm=100;
        tril_idx=find(tril(ones(size(stc,3)),-1));
        for perm_idx=1:n_perm
            for v_idx=1:size(stc,1)
                if(mod(v_idx,100)==0)
                    fprintf('permutation [%03d|%03d]...[%1.1f%%]...\r',perm_idx,n_perm,v_idx/size(stc,1)*100);
                end;
                
                
                for sequence_idx=1:length(sequence)
                    n_subj=size(data,2);
                    perm_swb_subj=randperm(n_subj);
                    perm_swb_subj=perm_swb_subj(1:round(length(perm_swb_subj)/2));
                    perm_swb_vec=ones(n_subj,1);
                    perm_swb_vec(perm_swb_subj)=-1;
                    perm_swb_mat=perm_swb_vec*perm_swb_vec';
                    
                    %data_perm=etc_phasescramble(data,'dim',1);
                    
                    tmp_perm=tril(squeeze(corr_buffer{hemi_idx}(sequence_idx,v_idx,:,:)).*perm_swb_mat,-1);
                    rr_perm=tmp_perm(tril_idx);
                    z_perm=0.5.*log((1+rr_perm)./(1-rr_perm))./(1/sqrt(size(data,1)/2.34-3));
                    %zm_perm=mean(z_perm);
                    %zmd_perm=median(z_perm);
                    
                    %cc(v_idx,:)=z(:)';
                    cc_perm{hemi_idx}(sequence_idx,:)=z_perm(:)';
                    %ccm_perm{hemi_idx}(sequence_idx)=zm_perm;
                    %ccmd_perm{hemi_idx}(sequence_idx)=zmd_perm;
                end;
                
                tmp=hypothesis{hypothesis_idx}(:)'*cc_perm{hemi_idx}(:,:);
                null_mean{hemi_idx}(v_idx,perm_idx)=mean(tmp);
                null_median{hemi_idx}(v_idx,perm_idx)=median(tmp);
            end;
        end;
        
        for v_idx=1:size(stc,1)
            zm_perm_p(v_idx)=length(find(null_mean{hemi_idx}(v_idx,:)>effect_mean{hemi_idx}(hypothesis_idx,v_idx)))./n_perm;
            zmd_perm_p(v_idx)=length(find(null_median{hemi_idx}(v_idx,:)>effect_median{hemi_idx}(hypothesis_idx,v_idx)))./n_perm;
        end;
        
        
        fprintf('\n');
        
        inverse_write_stc(repmat(zm_perm_p(:),[1 5]),a,b,c,sprintf('%s_gmr_mean_p-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(zmd_perm_p(:),[1 5]),a,b,c,sprintf('%s_gmr_median_p-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(1-zm_perm_p(:),[1 5]),a,b,c,sprintf('%s_gmr_mean_1mp-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(1-zmd_perm_p(:),[1 5]),a,b,c,sprintf('%s_gmr_median_1mp-%s.stc',oo,hemi{hemi_idx}));
    end;
end;
