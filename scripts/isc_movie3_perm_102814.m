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
root_path='/space_lin2/fhlin/chaplin/epi_data/chaplin_epi';
%root_path='/Users/fhlin_admin/workspace/chaplin/epi_data/chaplin_epi';

n_perm=1000;

perm_name='isc_102814_movie3_perm_z.mat';

output_stem='isc_102814_perm_z';

confound_polynomial_order=2;
confound_sinusoidal_order=3;

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
        
        for v_idx=1:size(stc,1)
            if(mod(v_idx,100)==0)
                fprintf('[%1.1f%%]...\r',v_idx/size(stc,1)*100);
            end;
            if(~isempty(D_prep))
                data=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:));
            else
                data=squeeze(stc(v_idx,:,:));
            end;
            tmp=tril(corrcoef(data),-1);
            rr=tmp(find(tmp));
            z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
            zm=median(z);
            cc(v_idx,:)=z(:)';
            
            ccm(v_idx)=zm;
        end;
        inverse_write_stc(cc,a,b,c,sprintf('%s-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(ccm(:),[1 5]),a,b,c,sprintf('%s_median-%s.stc',oo,hemi{hemi_idx}));
        
        data0=data;
        for perm_idx=1:n_perm
            fprintf('*');
            for v_idx=1:size(stc,1)
                %if(mod(perm_idx,10)==0)
                %    fprintf('permutation [%1.1f%%]...\r',perm_idx./n_perm.*100);
                %end;
                for ii=1:size(data,2)
                    data(:,ii)=surr_ft_algorithm3(data0(:,ii));
                end;
                tmp=tril(corrcoef(data),-1);
                rr=tmp(find(tmp));
                z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
                zm=median(z);
                cc(v_idx,:)=z(:)';
                
                
                if(v_idx==1&&perm_idx==1&&hemi_idx==1&&stc_idx==1)
                    effect_perm=zeros(length(zm),n_perm,2,length(stc_folder));
                end;
                effect_perm(v_idx,perm_idx,hemi_idx,stc_idx)=zm;
            end;
            
        end;
        fprintf('\n');
    end;
end;

save(perm_name,'effect_perm','-v7.3');
