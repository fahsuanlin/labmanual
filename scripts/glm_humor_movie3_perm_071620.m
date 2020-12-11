close all; clear all;

subject={
    '100701s1_FangWY';
    '100701s2_HsiaCC';
    '100706s1_ChengFH';
    '100706s2_PanYS';
    '100706s3_ChiJJ';
    '100708s1_HuangHJ';
    '100708s2_HuYS';
    '100708s3_LiuPY';
    '100713s1_KungYC';
    '100713s2_YangYC';
    '100713s3_WeiCC';
    'ChenHJ';
    'ChenTT';
    'ChuHF';
    'HuangPC';
    'LeeSM';
    'LinPF';
    'LiuJW';
    'ShenMH';
    'WuSY';
    };

stc_folder={
    'movie3a';
    'movie3b';
    %'movie5a';
    %'movie5b';
    %'movie7a';
    %'movie7b';
    };

hemi={
    'lh';
    'rh';
    };

%root_path='/space/maki5/1/users/fhlin/chaplin_fmri/epi_data';
root_path='/Users/fhlin_admin/workspace/chaplin/epi_data/chaplin_epi';

output_stem='glm_humor_071620';

confound_polynomial_order=2;
confound_sinusoidal_order=3;

n_perm=1000;
perm_name='glm_humor_071620_movie3_perm.mat';

file_score3='./interp_score_movie3_100914.mat';
file_score5='./interp_score_movie5_100914.mat';
file_score7='./interp_score_movie7_100914.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(file_score3);
for subj_idx=1:length(movie3)
    score{1}(subj_idx,:)=movie3(subj_idx).score_fmri;
end;
load(file_score5);
for subj_idx=1:length(movie5)
    score{2}(subj_idx,:)=movie5(subj_idx).score_fmri;
end;
load(file_score7);
for subj_idx=1:length(movie7)
    score{3}(subj_idx,:)=movie7(subj_idx).score_fmri;
end;


for stc_idx=1:length(stc_folder)
    oo=sprintf('%s_%s',output_stem,stc_folder{stc_idx});
    for hemi_idx=1:length(hemi)
        stc=[];
        for subj_idx=1:length(subject)
            d=dir(sprintf('%s/%s/%s/*_2_fsaverage*-%s.stc',root_path,subject{subj_idx},stc_folder{stc_idx},hemi{hemi_idx}));
            
            [stc(:,:,subj_idx),a,b,c]=inverse_read_stc(sprintf('%s/%s/%s/%s',root_path,subject{subj_idx},stc_folder{stc_idx},d(1).name));
            fprintf('[%s]::<%s>::<%s>\t%05d time ponits\n',stc_folder{stc_idx},subject{subj_idx},hemi{hemi_idx},size(stc,2));
            %			length(find(isnan(stc(:))))
        end;
        
        fprintf('\n');
        
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
        
        for subj_idx=1:size(stc,3)
            fprintf('subject [%d]...',subj_idx);
            data=squeeze(stc(:,:,subj_idx))'-D_prep*squeeze(stc(:,:,subj_idx))';
            contrast=squeeze(score{1}(subj_idx,:))';
            
            hrf=fmri_hdr([0:2:30]);
            hrf=hrf(:,1);
            contrast=(conv(contrast,hrf,'same'));
            
            beta=inv(contrast'*contrast)*contrast'*data;
            error=(data-contrast*beta);
            error_sig2=sum(error.^2,1)./(size(data,1)-size(contrast,2));
            df=size(data,1)-rank(contrast);
            
            cc=[1];
            t_stat(:,subj_idx)=cc'*beta./sqrt(error_sig2.*(cc'*inv(contrast'*contrast)*cc));
            effect(:,subj_idx)=cc'*beta;
            
            for perm_idx=1:n_perm
                %if(mod(perm_idx,10)==0)
                %    fprintf('permutation [%1.1f%%]...\r',perm_idx./n_perm.*100);
                %end;
                contrast=surr_ft_algorithm3(squeeze(score{1}(subj_idx,:)))';
                
                hrf=fmri_hdr([0:2:30]);
                hrf=hrf(:,1);
                contrast=(conv(contrast,hrf,'same'));
                
                beta=inv(contrast'*contrast)*contrast'*data;
                error=(data-contrast*beta);
                error_sig2=sum(error.^2,1)./(size(data,1)-size(contrast,2));
                df=size(data,1)-rank(contrast);
                
                cc=[1];
                %t_stat(:,subj_idx)=cc'*beta./sqrt(error_sig2.*(cc'*inv(contrast'*contrast)*cc));
                tmp=cc'*beta;
                if(subj_idx==1&&perm_idx==1&&hemi_idx==1&&stc_idx==1)
                    effect_perm=zeros(length(tmp),size(stc,3),n_perm,2,length(stc_folder));
                end;
                effect_perm(:,subj_idx,perm_idx,hemi_idx,stc_idx)=tmp;
            end;
            fprintf('\n');
        end;
        effect_z=mean(effect,2)./std(effect,0,2).*sqrt(size(stc,3));
        inverse_write_stc(repmat(effect_z(:),[1 5]),a,b,c,sprintf('%s_z-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(effect,a,b,c,sprintf('%s_effect-%s.stc',oo,hemi{hemi_idx}));
        

    end;
end;
save(perm_name,'effect_perm','-v7.3');
