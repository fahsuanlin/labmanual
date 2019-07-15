close all; clear all

%erf
file_erf='erf_041019.mat';

%fwd
file_fwd='meg_fwd_032819.mat';

output_stem='meg_041019';

SNR=100;
%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);
A2D=[];
for hemi_idx=1:2
    n_chan=size(A(hemi_idx).A,1);
    n_dip(hemi_idx)=size(A(hemi_idx).A,2);
    fprintf('[%d] MEG channels and [%d] dipoles\n',n_chan,n_dip(hemi_idx));
    n_source(hemi_idx)=n_dip(hemi_idx)/3;
    if(mod(n_dip(hemi_idx),3)~=0)
        
        fprintf('\n\n*** WARNING: The # of source is not 3-multiple! ****\n\n');
    else
        fprintf('[%d] sources\n',n_source(hemi_idx));
    end;
    Aa=reshape(A(hemi_idx).A,[n_chan, 3, n_source(hemi_idx)]);
    A_2d{hemi_idx}=reshape(Aa,[n_chan n_dip(hemi_idx)]);
    A2D=cat(2,A2D,A_2d{hemi_idx});
end;

R2D=[];
for hemi_idx=1:2
    R=ones(n_source(hemi_idx),1);
    %W=RA'*inv(A*R*A'+lambda.*C);
    RR=repmat(R,[1, 3])';
    R2D=cat(1,R2D,RR(:));
end;

RAt=repmat(R2D,[1 n_chan]).*(A2D');

ARAt=A2D*RAt;

p_signal=sum(diag(ARAt));

%noise covariance matrix
C=C(1).C;
p_noise=sum(diag(C));

lambda=p_signal/p_noise/SNR;

W2D=RAt*inv(ARAt+lambda.*C);

for trig_idx=1:size(erf_all,2)
    Y=erf_all(trig_idx).erf;
    baseline=find(erf_all(trig_idx).timeVec<0);
    Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
    
    Y_pred=zeros(size(Y));

    X_mne0=W2D*Y;
    X_mne=reshape(X_mne0,[3 sum(n_source) size(Y,2)]);
    X_mne=squeeze(sqrt(sum(X_mne.^2,1)));
    
    %Y_pred=A2D*X_mne0;
    
    X_dspm=bsxfun(@minus,X_mne, mean(X_mne(:,baseline),2));
    X_dspm=bsxfun(@rdivide,X_dspm,std(X_mne(:,baseline),0,2));
    
    t0=mean(diff(erf_all(trig_idx).timeVec)).*1e3; %sampling interval; ms;
    for hemi_idx=1:2
        switch hemi_idx
            case 1
                fn=sprintf('%s_%s-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(1:n_source(1),:);
            case 2
                fn=sprintf('%s_%s-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(n_source(1)+1:end,:);
        end;
        fprintf('\tsaving [%s]...\n',fn);
        inverse_write_stc(X_hemi,A(hemi_idx).v_idx,min(erf_all(trig_idx).timeVec).*1e3,t0,fn);
    end;
    
    %timeVec=erf_all(1).timeVec*1e3; %ms
    
    %load('meg_fwd_prep_032819.mat');
    %label=points_label(1:157);
    
    %h=etc_plotEF_kit(Y,'timeVec',timeVec,'label',label,'data_aux',Y_pred);
    
end;
