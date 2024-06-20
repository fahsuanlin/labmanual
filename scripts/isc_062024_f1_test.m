close all; clear all;

subject={
    'subj_01';
    'subj_02';
    'subj_03';
    'subj_04';
    'subj_05';
    'subj_06';
%     'subj_07';
%     'subj_08';
%     'subj_09';
%     'subj_10';
%     'subj_11';
%     'subj_12';
%     'subj_13';
%     'subj_14';
%     'subj_15';
%     'subj_16';
%     'subj_17';
%     'subj_18';
%     'subj_19';
    };


stc_folder={
    'f1';
    };

hemi={
    'lh';
    'rh';
    };

root_path='/space/maki6/1/users/fhlin/perspective';
root_path='/Users/fhlin/workspace/perspective';


output_stem='isc_062024_f1_z';

confound_polynomial_order=2;
confound_sinusoidal_order=3;

for stc_idx=1:length(stc_folder)
    oo=sprintf('%s',output_stem);
    for hemi_idx=1:length(hemi)
        stc=[];
        for subj_idx=1:length(subject)
            %d=dir(sprintf('%s/%s/%s/epi_data/unpack/*_2_fsaverage*-%s.stc',root_path,subject{subj_idx},stc_folder{stc_idx},hemi{hemi_idx}));
            fn=sprintf('%s_2_fsaverage_sfmcstc-%s.stc',subject{subj_idx},hemi{hemi_idx});
            [tmp,a,b,c]=inverse_read_stc(sprintf('%s/%s/epi_data/unpack/%s/%s',root_path,subject{subj_idx},stc_folder{stc_idx},fn));
            ll=size(tmp,2);
            tmp=tmp(:,4:ll-4);
            if(subj_idx>1)
                if(size(tmp,2)>size(stc,2))
                    tmp=tmp(:,1:end-1);
                end;
            end;
            stc(:,:,subj_idx)=tmp;
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
            %tmp=tril(corrcoef(squeeze(stc(v_idx,:,:))),-1);
            %cc(v_idx)=mean(tmp(find(tmp)));
            %for subj_idx=1:size(stc,3)
            if(~isempty(D_prep))
                data=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:));
            else
                data=squeeze(stc(v_idx,:,:));
            end;
            %end;
            %data=detrend(squeeze(stc(v_idx,:,:)));
            %[u,s,v]=svd(data);
            %s=diag(s);
            %cc(v_idx)=s(1).^2./sum(s.^2);
            %             for s_idx=1:size(data,2)
            %                 surr_idx=setdiff([1:size(data,2)],s_idx);
            %                 t1=data(:,s_idx);
            %                 t2=mean(data(:,surr_idx),2);
            %                 tmp_corr=corrcoef(t1,t2);
            %                 rr(s_idx)=tmp_corr(2,1);
            %             end;
            tmp=tril(corrcoef(data),-1);
            rr=tmp(find(tmp));
            z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
            zm=median(z);
            %chi2=sum(log(1-normcdf(abs(z)))).*(-2); %http://en.wikipedia.org/wiki/Fisher%27s_method
            %chi2m=sum(log(1-normcdf(abs(zm)))).*(-2); %http://en.wikipedia.org/wiki/Fisher%27s_method
            %if(~isempty(find(chi2(:)<eps)))
            %	fprintf('Chi2 = 0!!\n');
            %	keyboard;
            %end;
            %if(~isempty(find(chi2(:)>1e3)))
            %        fprintf('Chi2 = Inf!!\n');
            %        keyboard;
            %end;
            
            %chi2(find(isnan(chi2)))=0;
            %chi2(find(isinf(chi2)))=1000;
            cc(v_idx,:)=z(:)';
            
            %chi2m(find(isnan(chi2m)))=0;
            %chi2m(find(isinf(chi2m)))=1000;
            ccm(v_idx)=zm;
        end;
        
        inverse_write_stc(cc,a,b,c,sprintf('%s-%s.stc',oo,hemi{hemi_idx}));
        inverse_write_stc(repmat(ccm(:),[1 5]),a,b,c,sprintf('%s_median-%s.stc',oo,hemi{hemi_idx}));
    end;
end;

