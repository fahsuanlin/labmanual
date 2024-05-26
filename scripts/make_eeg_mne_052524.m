close all; clear all;

file_fwd='eeg_fwd_052524.mat';
file_bem='fwd_prep_041618.mat'; %including the source definition

freqVec=[15];

flag_morph=0;
subject='180411_PYY';
target_subject='fsaverage';

output_stem={
    'SSVEP_out_1_042124';
    };

file_erp={
    'erp_avg_outside.mat';
    };

load(file_bem);

for f_idx=1:length(file_erp)
    %prepare ERP
    load(file_erp{f_idx});


    bad_channel=[];

    Y=erp_avg{3}.erp; %<----the 3rd condition in the erp_avg variable.
    timeVec=erp_avg{3}.timeVec; 
    label=erp_avg{3}.electrode_name;
    %noise covariance matrix


    %calculate inverse
    load(file_fwd);
    F=cat(2,A(1).A,A(2).A);


    %find correspondence channels
    [d,i1,i2]=intersect(A(1).label, lower(label));
    F=F(i1,:);
    Y=Y(i2,:);

    F(bad_channel,:)=[];
    Y(bad_channel,:)=[];


    %noise covariance matrix
    baseline_idx=[1:round(length(timeVec)./1e3):length(timeVec)];
    C=cov(Y(:,baseline_idx).'); %noise covariance


    W=F'*inv(F*F'+trace(F*F')/trace(C)*0.1.*C); %MNE


    X=W*Y;
    X_mne=reshape(X,[3,size(X,1)./3,size(X,2)]);
    X_mne_abs=squeeze(sqrt(sum(X_mne.^2,1)));

    inverse_write_stc(X_mne_abs(1:src(1).nuse,:),src(1).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_mne_abs-lh.stc',output_stem{f_idx}));
    inverse_write_stc(X_mne_abs(src(1).nuse+1:end,:),src(2).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_mne_abs-rh.stc',output_stem{f_idx}));

    baseline_idx=find(timeVec<0);
    X_dspm=bsxfun(@minus,X_mne_abs,mean(X_mne_abs(:,baseline_idx),2));
    X_dspm=bsxfun(@rdivide,X_dspm,std(X_mne_abs(:,baseline_idx),0,2));

    inverse_write_stc(X_dspm(1:src(1).nuse,:),src(1).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_dspm-lh.stc',output_stem{f_idx}));
    inverse_write_stc(X_dspm(src(1).nuse+1:end,:),src(2).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_dspm-rh.stc',output_stem{f_idx}));
       


    if(~isempty(freqVec))
        for freq_idx=1:length(freqVec)
            tfr=abs(inverse_waveletcoef(freqVec(freq_idx),X,1/mean(diff(timeVec)),5));

            tfr_mne=reshape(tfr,[3,size(tfr,1)./3,size(tfr,2)]);
            tfr_mne_abs=squeeze(sqrt(sum(tfr_mne.^2,1)));
   
            tfr_dspm=bsxfun(@minus,tfr_mne_abs,mean(tfr_mne_abs(:,baseline_idx),2));
            tfr_dspm=bsxfun(@rdivide,tfr_dspm,std(tfr_mne_abs(:,baseline_idx),0,2));
            
            inverse_write_stc(tfr_dspm(1:src(1).nuse,:),src(1).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_%2.1fhz_dspm-lh.stc',output_stem{f_idx},freqVec(freq_idx)));
            inverse_write_stc(tfr_dspm(src(1).nuse+1:end,:),src(2).vertno,min(timeVec).*1e3,mean(diff(timeVec)).*1e3,sprintf('%s_%2.1fhz_dspm-rh.stc',output_stem{f_idx},freqVec(freq_idx)));
        end;
    end;
end;
