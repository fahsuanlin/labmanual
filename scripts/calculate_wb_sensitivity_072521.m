close all; clear all

%ncov
file_ncov='erf_050719.mat';

%fwd
file_fwd='seeg_fwd_wb_dec_072521.mat';

output_stem='seeg_wb_sensitivity_072521';

SNR=10;
%%%%%%%%%%%%%%

load(file_ncov);

load(file_fwd);

%search electrode correspondence.
%only electrodes in the noise covariance matrix are included!
%the variable 'A_idx(j)' is the index of electrodes in the forward solution
%that matches the electrodes j in the noise covariance matrix
%
c_keep=ones(1, length(C(1).name));

for d_idx=1:length(C(1).name)
    n=C(1).name{d_idx};
    n=n(2:end);
    n=erase(n,"'");
    
    if(sum(abs(strcmp(A(1).name, n)))>eps)
        %find the matched lead fields in the forward solution
        A_idx(d_idx)=find(strcmp(A(1).name, n));
    else
        fprintf('contact [%s] in the noise covariance matrix not found in the forward solution!\n',n);
        c_keep(d_idx)=0;
    end;
end;


%noise covariance matrix
for idx=1:length(C)
    if(idx==1) Ctmp=C(idx).C; else Ctmp=Ctmp+C(idx).C; end;
end;
C=Ctmp./length(C);


C=C(find(c_keep),find(c_keep));
%C=diag(diag(C));

[uc,sc]=svd(C);
ss=diag(sc);
idx=find(cumsum(ss.^2)./sum(ss.^2)<0.99);
%idx=find(diag(sc)>eps);
tmp=diag(sc);
sc_half=zeros(1,size(C,1));
sc_half(idx)=1./sqrt(tmp(idx));
sc_half=diag(sc_half);
C_half=uc*sc_half;


A2D=[];
A2D_sens=[];
A2D_sens_nn=[];
for hemi_idx=1:2
    %A_idx=[1:size(A(hemi_idx).A,1)];
    if(isfield(A(hemi_idx),'wb_loc'))
        %whole-brain forward solution
        source_idx=[1:size(A(hemi_idx).A,2)];
    else
        %cortical source space forward solution        
        source_idx=[1:size(A(hemi_idx).loc,1)];
    end;

    %make sure the lead field matched the data
    %A(hemi_idx).A=A(hemi_idx).A(A_idx,source_idx); %now the order of the electrode is matched to that in the noise covariance matrix
    A(hemi_idx).A=A(hemi_idx).A(A_idx,:); %now the order of the electrode is matched to that in the noise covariance matrix
    
    %A(hemi_idx).A=A(hemi_idx).A(:,1:length(A(hemi_idx).v_idx)*3);
    
    n_chan=size(A(hemi_idx).A,1);
    n_dip(hemi_idx)=size(A(hemi_idx).A,2);
    fprintf('[%d] SEEG contacts and [%d] dipoles\n',n_chan,n_dip(hemi_idx));
    n_source(hemi_idx)=n_dip(hemi_idx)/3;
    if(mod(n_dip(hemi_idx),3)~=0)
        
        fprintf('\n\n*** WARNING: The # of source is not 3-multiple! ****\n\n');
    else
        fprintf('[%d] sources\n',n_source(hemi_idx));
    end;
    Aa=reshape(A(hemi_idx).A,[n_chan, 3, n_source(hemi_idx)]);
    A_2d{hemi_idx}=reshape(Aa,[n_chan n_dip(hemi_idx)]);
    A2D=cat(2,A2D,A_2d{hemi_idx});
    A2D_sens=cat(1,A2D_sens,squeeze(sum(sum(Aa.^2,1),2)));
    
    tmp=squeeze(sum((C_half'*A_2d{hemi_idx}).^2,1));
    tmp=reshape(tmp,[3 n_source(hemi_idx)]);
    tmp=sum(tmp,1);
    A2D_sens_nn=cat(1,A2D_sens_nn,tmp(:));
end;
%A2D_sens_nn=A2D_sens_nn./max(A2D_sens_nn).*1;

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
p_noise=sum(diag(C));

lambda=p_signal/p_noise/SNR;

W2D=RAt*inv(ARAt+lambda.*C);

%for trig_idx=1:size(erf_all,2)
for trig_idx=1:1
    %Y=erf_all(trig_idx).erf;
    %baseline=find(erf_all(trig_idx).timeVec<0);
    %Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
    
    %noise_power=mean(sum(Y(:,baseline).^2,1));
    noise_power=sum(diag(C));
    signal_power=mean(A2D_sens);
    %Y=Y./sqrt(noise_power./signal_power);
    
    %X_mne0=W2D*Y(:,baseline);
   
    %X_mne=reshape(X_mne0,[3 sum(n_source) size(Y(:,baseline),2)]);
    %X_mne=squeeze(sqrt(sum(X_mne.^2,1)));
    %clear X_mne0;
    %nn=sum(X_mne,2)./size(X_mne,2);
    %SNR_map=10*log(A2D_sens./nn); %in dB
    S_map=A2D_sens; 
    SNR_map=A2D_sens_nn;
    
    %t0=mean(diff(erf_all(trig_idx).timeVec)); %sampling interval; ms;
    for hemi_idx=1:2
        switch hemi_idx
            case 1
                fn=sprintf('%s_snr-lh.stc',output_stem);
                X_hemi=repmat(SNR_map(1:length(A(hemi_idx).v_idx)),[1 5]);
                
            case 2
                fn=sprintf('%s_snr-rh.stc',output_stem);
                X_hemi=repmat(SNR_map(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx)),[1 5]);
        end;
        fprintf('\tsaving [%s]...\n',fn);
        inverse_write_stc(X_hemi,A(hemi_idx).v_idx,0,1,fn);

        switch hemi_idx
            case 1
                fn=sprintf('%s_s-lh.stc',output_stem);
                X_hemi=repmat(S_map(1:length(A(hemi_idx).v_idx)),[1 5]);
                
            case 2
                fn=sprintf('%s_s-rh.stc',output_stem);
                X_hemi=repmat(S_map(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx)),[1 5]);
        end;
        fprintf('\tsaving [%s]...\n',fn);
        inverse_write_stc(X_hemi,A(hemi_idx).v_idx,0,1,fn);
    end;
    
    fn=sprintf('%s_snr-vol.stc',output_stem);
    fprintf('\tsaving [%s]...\n',fn);
    inverse_write_stc(repmat(SNR_map(:),[1 5]),[0:size(SNR_map,1)-1],0,1,fn);    


    fn=sprintf('%s_s-vol.stc',output_stem);
    fprintf('\tsaving [%s]...\n',fn);
    inverse_write_stc(repmat(S_map(:),[1 5]),[0:size(S_map,1)-1],0,1,fn);    
end;
