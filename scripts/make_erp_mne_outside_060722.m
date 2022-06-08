close all; clear all;

file_fwd='eeg_fwd_wb_060722.mat';
file_bem='eeg_fwd_prep_060722.mat'; %including the source definition

bad_channel=[32 33 34]; %    'Ref', 'Gnd', 'REF.'

flag_morph=0;
subject='180322_PYW';
target_subject='fsaverage';

sf_target=500; %target sampling rate; Hz

output_stem='erp_mne_outside_060722';

file_erp={
    'erp_outside.mat';
    };
erp_trigger={[3]};

SNR=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_bem);
load(file_fwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for f_idx=1:length(file_erp)
    %prepare ERP
    load(file_erp{f_idx});
    
    tmp=[];
    
    for trigger_idx=1:length(erp_trigger)
        erp_counter=0;
        
        for i=1:size(erp,1)
            if(isempty(tmp))
                tmp=erp{i,erp_trigger{trigger_idx}}.erp;
            else
                tmp=tmp+erp{i,erp_trigger{trigger_idx}}.erp;
            end;
            erp_counter=erp_counter+1;
        end;
        
        
        Y=tmp./erp_counter;
        timeVec=erp{erp_trigger{trigger_idx}}.timeVec;
        
        
        %down-sampling
        sf0=1./mean(diff(timeVec));
        oversampling_ratio=ceil(sf0./sf_target);
        for ch_idx=1:size(Y,1)
            Y(ch_idx,:)=filtfilt(ones(oversampling_ratio,1),1,Y(ch_idx,:));
        end;
        Y=Y(:,1:oversampling_ratio:end);
        timeVec=timeVec(1:oversampling_ratio:end);
        
        
        %noise covariance matrix; based on evoked response
        baseline_idx=find(timeVec<0);
        C=cov(Y(:,baseline_idx).'); %noise covariance
        
        data_idx_remove=[];
        %search correspondence
        A_idx_tmp=zeros(1, length(erp{1}.electrode_name)).*nan;
        for d_idx=1:length(erp{1}.electrode_name)
            n=erp{1}.electrode_name{d_idx};
            
            
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
            RR=repmat(R,[1, 3])';
            R2D=cat(1,R2D,RR(:));
        end;
        
        RAt=repmat(R2D,[1 n_chan]).*(A2D');
        
        ARAt=A2D*RAt;
        
        p_signal=sum(diag(ARAt));
        
        C=diag(diag(C));
        p_noise=sum(diag(C));
        
        lambda=p_signal/p_noise/SNR;
        
        W2D=RAt*inv(ARAt+lambda.*C);
        
        baseline=find(timeVec<0);
        Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
        
        X_mne=W2D*Y;
        X_mne3=reshape(X_mne,[3,size(X_mne,1)/3,size(X_mne,2)]);
        
        X_dspm=etc_z(X_mne,find(timeVec<0));
        
        t0=mean(diff(timeVec)); %sampling interval; ms;
        
        for dd=1:3 %x,y,z
            switch dd
                case 1
                    dd_str='x';
                case 2
                    dd_str='y';
                case 3
                    dd_str='z';
            end;
            
            for hemi_idx=1:2
                switch hemi_idx
                    case 1
                        fn_mne=sprintf('%s_%03d_mne_%s-lh.stc',output_stem,erp_trigger{trigger_idx},dd_str);
                        X_hemi_mne=X_mne(dd:3:length(A(hemi_idx).v_idx)*3,:);
                        fn_dspm=sprintf('%s_%03d_dspm_%s-lh.stc',output_stem,erp_trigger{trigger_idx},dd_str);
                        X_hemi_dspm=X_dspm(dd:3:length(A(hemi_idx).v_idx)*3,:);
                        hemi='lh';
                    case 2
                        fn_mne=sprintf('%s_%03d_mne_%s-rh.stc',output_stem,erp_trigger{trigger_idx},dd_str);
                        X_hemi_mne=X_mne(n_source(1)*3+dd:3:n_source(1)*3+length(A(hemi_idx).v_idx)*3,:);
                        fn_dspm=sprintf('%s_%03d_dspm_%s-rh.stc',output_stem,erp_trigger{trigger_idx},dd_str);
                        X_hemi_dspm=X_dspm(n_source(1)*3+dd:3:n_source(1)*3+length(A(hemi_idx).v_idx)*3,:);
                        hemi='rh';
                end;
                fprintf('\tsaving [%s]...\n',fn_mne);
                inverse_write_stc(X_hemi_mne,A(hemi_idx).v_idx,min(timeVec),t0,fn_mne);
                inverse_write_stc(X_hemi_dspm,A(hemi_idx).v_idx,min(timeVec),t0,fn_dspm);
            end;
            
            
        end;
        
        
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    fn_mne=sprintf('%s_%03d_mne_%s-lh.stc',output_stem,erp_trigger{trigger_idx},'cc');
                    X_hemi_mne=squeeze(sum(repmat(A(hemi_idx).ori',[1 1 size(X_mne3,3)]).*X_mne3(:,1:length(A(1).v_idx),:),1));
                    
                    fn_dspm=sprintf('%s_%03d_dspm_%s-lh.stc',output_stem,erp_trigger{trigger_idx},'cc');
                    X_hemi_dspm=etc_z(X_hemi_mne,find(timeVec<0));
                    hemi='lh';
                case 2
                    fn_mne=sprintf('%s_%03d_mne_%s-rh.stc',output_stem,erp_trigger{trigger_idx},'cc');
                    X_hemi_mne=squeeze(sum(repmat(A(hemi_idx).ori',[1 1 size(X_mne3,3)]).*X_mne3(:,size(A(1).A,2)./3+1:size(A(1).A,2)./3+length(A(2).v_idx),:),1));
                    
                    fn_dspm=sprintf('%s_%03d_dspm_%s-rh.stc',output_stem,erp_trigger{trigger_idx},'cc');
                    X_hemi_dspm=etc_z(X_hemi_mne,find(timeVec<0));
                    hemi='rh';
            end;
            fprintf('\tsaving [%s]...\n',fn_mne);
            inverse_write_stc(X_hemi_mne,A(hemi_idx).v_idx,min(timeVec),t0,fn_mne);
            inverse_write_stc(X_hemi_dspm,A(hemi_idx).v_idx,min(timeVec),t0,fn_dspm);
        end;
        
        %morphing STC....
        for dd=1:3 %x,y,z
            switch dd
                case 1
                    dd_str='x';
                case 2
                    dd_str='y';
                case 3
                    dd_str='z';
            end;
            for hemi_idx=1:2
                switch hemi_idx
                    case 1
                        hemi='lh';
                    case 2
                        hemi='rh';
                end;
                fn_in=sprintf('%s_%03d_mne_%s-%s.stc',output_stem,erp_trigger{trigger_idx},dd_str,hemi);
                fn_out=sprintf('%s_2_%s_%03d_%s_mne_%s',subject,target_subject,output_stem,erp_trigger{trigger_idx},dd_str);
                cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_in, target_subject, fn_out, hemi);
                eval(cmd);
            end;
        end;
    end;
end;
return;

