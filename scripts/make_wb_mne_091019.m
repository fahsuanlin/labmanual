close all; clear all

%erf
file_erf='erf_050719.mat';
flag_auto_remove_comma=1; %remove ' from the electrode name automatically

%fwd
file_fwd='seeg_fwd_wb_091019.mat';

output_stem='seeg_wb_mne_091019';

subject='s031';
target_subject='fsaverage';

SNR=100;

flag_morph=0;
%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);

data_idx_remove=[];
%search correspondence
A_idx_tmp=zeros(1, length(erf_all(1).name)).*nan;
for d_idx=1:length(erf_all(1).name)
    n=erf_all(1).name{d_idx};
    n=n(2:end);
    if(flag_auto_remove_comma)
        n=erase(n,"'");
    end;
    
    %find the matched lead fields in the forward solution
    ii=find(strcmp(A(1).name, n));
    if(~isempty(ii))
        A_idx_tmp(d_idx)=find(strcmp(A(1).name, n));
    else
        data_idx_remove(end+1)=d_idx;
        fprintf('channel [%s] in data has no corresponding entry in the electrode list!\n',n);
    end;
end;
ii=(find(~isnan(A_idx_tmp)));
A_idx=A_idx_tmp(ii);

A2D=[];
for hemi_idx=1:2
    %make sure the lead field matched the data
    A(hemi_idx).A=A(hemi_idx).A(A_idx,:);
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
for idx=1:length(erf_all)
    if(idx==1) Ctmp=C(idx).C; else Ctmp=Ctmp+C(idx).C; end;
end;
C=Ctmp./length(erf_all);
C(data_idx_remove,:)=[];
C(:,data_idx_remove)=[];

C=diag(diag(C));
p_noise=sum(diag(C));

lambda=p_signal/p_noise/SNR;

W2D=RAt*inv(ARAt+lambda.*C);

for trig_idx=1:size(erf_all,2)
    Y=erf_all(trig_idx).erf;
    baseline=find(erf_all(trig_idx).timeVec<0);
    Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
    
    Y_pred=zeros(size(Y));
    
    Y(data_idx_remove,:)=[];
    
    X_mne0=W2D*Y;
    
    X_mne=reshape(X_mne0,[3 sum(n_source) size(Y,2)]);
    X_mne=squeeze(sqrt(sum(X_mne.^2,1)));
    clear X_mne0;
    
    X_dspm=bsxfun(@minus,X_mne, mean(X_mne(:,baseline),2));
    X_dspm=bsxfun(@rdivide,X_dspm,std(X_dspm(:,baseline),0,2));
    
    t0=mean(diff(erf_all(trig_idx).timeVec)); %sampling interval; ms;
    for hemi_idx=1:2
        switch hemi_idx
            case 1
                hemi='lh';
                fn=sprintf('%s_%s-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(1:length(A(hemi_idx).v_idx),:);
                fn_mne=sprintf('%s_%s_mne-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi_mne=X_mne(1:length(A(hemi_idx).v_idx),:);
            case 2
                hemi='rh';
                fn=sprintf('%s_%s-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx),:);
                fn_mne=sprintf('%s_%s_mne-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi_mne=X_mne(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx),:);
        end;
        fprintf('\tsaving [%s]...\n',fn);
        inverse_write_stc(X_hemi,A(hemi_idx).v_idx,min(erf_all(trig_idx).timeVec),t0,fn);
        fprintf('\tsaving [%s]...\n',fn_mne);
        inverse_write_stc(X_hemi_mne,A(hemi_idx).v_idx,min(erf_all(trig_idx).timeVec),t0,fn_mne);
        
    end;

    if(flag_morph)
        %morphing...
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    hemi='lh';
                case 2
                    hemi='rh';
            end;
            fn_out=sprintf('%s_2_%s_%s_%s-%s.stc',subject,target_subject,output_stem,erf_all(trig_idx).trig_str,hemi);
            cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn, target_subject, fn_out, hemi);
            eval(cmd);
            fn_mne_out=sprintf('%s_2_%s_%s_%s_mne-%s.stc',subject,target_subject,output_stem,erf_all(trig_idx).trig_str,hemi);
            cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_mne, target_subject, fn_mne_out, hemi);
            eval(cmd);
        end;
    end;
    fn=sprintf('%s_%s-vol.stc',output_stem,erf_all(trig_idx).trig_str);
    fprintf('\tsaving [%s]...\n',fn);
    inverse_write_stc(X_dspm,[0:size(X_dspm,1)-1],min(erf_all(trig_idx).timeVec),t0,fn);
    
    fn=sprintf('%s_%s_mne-vol.stc',output_stem,erf_all(trig_idx).trig_str);
    fprintf('\tsaving [%s]...\n',fn);
    inverse_write_stc(X_mne,[0:size(X_mne,1)-1],min(erf_all(trig_idx).timeVec),t0,fn);
    
    
end;
