close all; clear all;

headerFile = {
    '../eeg_raw/SSVEP_noMR_1.vhdr';
    };

markerFile={
    '../eeg_raw/SSVEP_noMR_1.vmrk';
    };

file_suffix='042124';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG setup
%
n_chan=32; %64-channel EEG;
select_channel={ 'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'    'ECG'};
eeg_channel=[1:31];
ecg_channel=[32];

TR=2; %second

%these two tokens are required (but values are arbitrary) for data
%collected inside MRI
trigger_token=1e3;
sync_token=1e2;

flag_reref= [0 0 0 0 0 0]; %referencing: remove the average time course across all electrodes from the time course at each individual electrode
flag_hp=    [1 1 1 1 1 1]; %high-pass filtering at 0.1 Hz
flag_lp=    [1 1 1 1 1 1]; %low-pass filtering at 70 Hz
flag_aas=   [0 0 0 0 0 0]; %AAS for gradient artifact suppression
flag_bcg=    [1 1 1 1 1 1]; %BCG artifact suppression
flag_bcg_ccm =  [0 0 0 0 0 0]; %BCG artiface suppression using DMH
time_trim=[];     %trim off the first few seconds...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f_idx=1:length(headerFile)
    [dummy,fstem]=fileparts(headerFile{f_idx});
    fstem=sprintf('%s_%s',fstem,file_suffix);
    fprintf('reading [%s]...\n',fstem);
    
    % first get the continuous data as a matlab array
    eeg{f_idx} = double(bva_loadeeg(headerFile{f_idx}));
    
    % meta information such as samplingRate (fs), labels, etc
    [fs(f_idx) label meta] = bva_readheader(headerFile{f_idx});
    
    found_channel={};
    for s_idx=1:length(select_channel)
        IndexC = strcmp(lower(label),lower(select_channel{s_idx})); %change all labels into lower case
        Index = find(IndexC);
        
        if(~isempty(Index))
            fprintf('\tChannel [%s] found:: index=%03d \r',select_channel{s_idx},Index);
            if(strcmp(lower(select_channel{s_idx}),'ecg'))
                ecg_channel=Index;
            else
                eeg_channel(s_idx)=Index;
            end;
            found_channel{end+1}=select_channel{s_idx};
        else
            fprintf('\tChannel [%s] not found! \r',select_channel{s_idx});
        end;
    end;
    fprintf('\n');
    
    ecg_orig=eeg{f_idx}(ecg_channel,:)';
    eeg_orig=eeg{f_idx}(eeg_channel,:);
    ecg{f_idx}=eeg{f_idx}(ecg_channel,:)';
    eeg{f_idx}=eeg{f_idx}(eeg_channel,:);
    
    
    %read maker file
    fprintf('\treading triggers...\n');
    trigger{f_idx}=etc_read_vmrk(markerFile{f_idx});
    
    
    EEG_before_aas=eeg{f_idx};
    
    if(flag_aas(f_idx))
        fprintf('\tAAS...\n');
        gradient_trigger{f_idx}=zeros(size(eeg{f_idx},2),1);
        iidx=find(trigger{f_idx}.event==trigger_token);
        gradient_trigger{f_idx}(trigger{f_idx}.time(iidx))=1e3;
        ecg{f_idx}=eeg_ga(ecg{f_idx},gradient_trigger{f_idx},TR,fs(f_idx),'flag_display',0,'flag_ma_aas',1,'flag_aas_svd',0,'flag_anchor_bnd',0,'n_ma_aas',7); %AAS on ECG
        ecg{f_idx}=sgolayfilt(ecg{f_idx},6,301); %smoothing out residual GA in ECG
        eeg{f_idx}=eeg_ga(eeg{f_idx},gradient_trigger{f_idx},TR,fs(f_idx),'flag_display',0,'flag_ma_aas',1,'flag_aas_svd',0,'flag_anchor_bnd',0,'n_ma_aas',7); %AAS on EEG
    end;
    
    
    
    %high-pass filtering
    if(flag_hp(f_idx))
        fprintf('\tHP...\n');
        nn=round(fs(f_idx).*1); %1 Hz
        for ch_idx=1:size(eeg{f_idx},1)
            fprintf('*');
            lp_data=filtfilt(ones(1,nn)./nn,1,eeg{f_idx}(ch_idx,:));
            eeg{f_idx}(ch_idx,:)=eeg{f_idx}(ch_idx,:)-lp_data;
        end;
%         if(~isempty(eeg_ref{f_idx}))
%             fprintf('.');
%             lp_data=filtfilt(ones(1,nn)./nn,1,eeg_ref{f_idx}(:)');
%             eeg_ref{f_idx}=eeg_ref{f_idx}(:)'-lp_data;
%         end;
        fprintf('~');
        lp_data=filtfilt(ones(1,nn)./nn,1,ecg{f_idx}(:)');
        ecg{f_idx}=ecg{f_idx}(:)'-lp_data;
        fprintf('\n');
%         
%         %high-pass filtering (0.1 Hz)
%         Wn = 0.1*2/fs(f_idx);
%         N = 3; % order of 3 less processing
%         [a,b] = butter(N,Wn,'high'); %bandpass filtering
%         for ch_idx=1:size(eeg{f_idx},1)
%             eeg{f_idx}(ch_idx,:) = filtfilt(a,b,eeg{f_idx}(ch_idx,:));
%         end;
%         if(~isempty(eeg_ref{f_idx}))
%             eeg_ref{f_idx}=filtfilt(a,b,eeg_ref{f_idx}(:))';
%         end;
%         
%         ecg{f_idx}=filtfilt(a,b,ecg{f_idx}(:))';
    else
        ecg{f_idx}=ecg{f_idx}(:)';
    end;      
    
    EEG_before_bcg=eeg{f_idx};
    
    %remove first few seconds (if needed....)
    if(~isempty(time_trim))
        time_trim_idx=round(time_trim*fs(f_idx));
        
        fprintf('trimming [%1.1f] s data {(%d) samples}....\n',time_trim, time_trim_idx);
        ecg{f_idx}=ecg{f_idx}(time_trim_idx+1:end);
        eeg{f_idx}=eeg{f_idx}(:,time_trim_idx+1:end);
        %trigger{f_idx}.time=trigger{f_idx}.time-time_trim_idx;
        
        %%update trigger info
        %idx=find(trigger{f_idx}.time<0);
        %trigger{f_idx}.time(idx)=[];
        %trigger{f_idx}.event(idx)=[];
    end;
    
    
   
    
    [qrs_amp_raw,qrs{f_idx}]=pan_tompkin2(ecg{f_idx},fs(f_idx));
    %%plot ECG and QRS detection
    %tt=[1:length(ecg{f_idx})]./fs(f_idx);
    %figure; plot(tt,ecg{f_idx}); hold on;
    %line(repmat(tt(qrs{f_idx}),[2 1]),repmat([min(ecg{f_idx}-mean(ecg{f_idx}))/2; max(ecg{f_idx}-mean(ecg{f_idx}))/2],size(qrs{f_idx})),'LineWidth',2.5,'LineStyle','-.','Color','r');
    
    eeg_bcg_pred{f_idx}=[];
    if(flag_bcg(f_idx))
        fprintf('\tBCG...\n');
        %[eeg{f_idx}, qrs_i_raw, bcg_all, ecg_all, bad_trials]=eeg_bcg(eeg{f_idx},ecg{f_idx},fs(f_idx),'flag_display',0,'bcg_nsvd',3,'n_ma_bcg',20,'flag_anchor_ends',0,'flag_post_ssp',0,'flag_badrejection',0); %BCG correction on EEG
        
        iqi=diff(qrs{f_idx})./fs(f_idx);
        iqi=sort(iqi);
        iqi_min=iqi(round(length(iqi).*0.5));
        tpre=0.2.*iqi_min;
        tpost=0.8.*iqi_min;
        
        %down-sampling
        tmp=[];
        for ch_idx=1:size(eeg{f_idx},1)
            tmp0=filtfilt(ones(100,1)./100,1,eeg{f_idx}(ch_idx,:));
            tmp(ch_idx,:)=tmp0(1:100:end);
        end;
        eeg_dec=tmp;
        
        tmp=filtfilt(ones(100,1)./100,1,ecg{f_idx});
        ecg_dec=tmp(1:100:end);
        
        fs_dec=fs(f_idx)/100;
        
        eeg_dec_bcg_pred=[];
        if(~flag_bcg_ccm(f_idx))
            [eeg_dec_bcg, qrs_i_raw, bcg_all, ecg_all, bad_trials]=eeg_bcg(eeg_dec,ecg_dec,fs_dec,'flag_display',0,'bcg_nsvd',3,'n_ma_bcg',21,'flag_anchor_ends',0,'flag_post_ssp',0,'flag_badrejection',0,'BCG_tPre',tpre,'BCG_tPost',tpost,'flag_pan_tompkin2',1); %BCG correction on EEG
        else
            e=8;
            tau=1;
            nn=10;
            [eeg_dec_bcg, qrs_i_raw, eeg_dec_bcg_pred]=eeg_bcg_ccm2(eeg_dec,ecg_dec,fs_dec,'e',e,'tau',tau,'nn',nn,'flag_display',0,'n_ecg',30);

            %[eeg_dec_bcg, eeg_dec_bcg_pred]=eeg_bcg_ccm5(eeg_dec,ecg_dec,fs_dec,'flag_display',0,'flag_auto_hp',1,'flag_pan_tompkin2',1);
        end;

        %up-samping
        for ch_idx=1:size(eeg_dec_bcg,1)
            tmp=interp(eeg_dec_bcg(ch_idx,:),100);
            %eeg{f_idx}(ch_idx,:)=sgolayfilt(tmp,4,301);
            eeg{f_idx}(ch_idx,:)=tmp(1:size(eeg{f_idx},2));
            
            if(~isempty(eeg_dec_bcg_pred))
                tmp=interp(eeg_dec_bcg_pred(ch_idx,:),100);
                eeg_bcg_pred{f_idx}(ch_idx,:)=tmp(1:size(eeg{f_idx},2));
            else
                eeg_bcg_pred{f_idx}=[];
            end;
        end;
        
    end;
    
    
    EEG_after_bcg=eeg{f_idx};

    %re-referencing
    if(flag_reref(f_idx))
        fprintf('\treferencing...\n');
        eeg_ref{f_idx}=mean(eeg{f_idx},1);
        for ch_idx=1:size(eeg{f_idx},1)
            eeg{f_idx}(ch_idx,:)=eeg{f_idx}(ch_idx,:)-eeg_ref{f_idx};
        end;
    else
        eeg_ref{f_idx}=[];
    end;

    EEG_after_reref=eeg{f_idx};

    
    %assign a trigger to R-peaks of ECG
    trigger_ecg{f_idx}.event=33.*ones(1,length(qrs{f_idx}(:)));
    trigger_ecg{f_idx}.time=qrs{f_idx}(:)';
    %if data has been trimmed, ECG trigger should be adjusted.
    if(~isempty(time_trim))
        trigger_ecg{f_idx}.time=trigger_ecg{f_idx}.time+round(fs(f_idx).*time_trim);
    end;
    %merge ECG triggers
    trigger{end}=etc_trigger_append(trigger{end},trigger_ecg{f_idx});
    
    
    
    %update trigger information by including ECG triggers
    TRIGGER_ECG=trigger_ecg{end};
    TRIGGER=trigger{end};
    
    EEG=eeg{end}; %<--------
    ECG=ecg{end};
    
    if(~isempty(eeg_bcg_pred{f_idx}))
        EEG_bcg_pred=eeg_bcg_pred{f_idx};
    else
        EEG_bcg_pred=[];
    end;
    
    %low-poass filtering EEG (cut-off = 70 Hz)
    if(flag_lp(f_idx))
        a=1;
        b=ones(round(fs(f_idx)./70),1)./round(fs(f_idx)./70);
        for idx=1:size(EEG,1)
            EEG(idx,:)=filtfilt(b,a,EEG(idx,:));
        end;
    end;
    
    %compensate the EEC/ECG data with zeros for the truncated time, if any
    EEG=cat(2,eeg_orig(:,1:fs(f_idx).*time_trim),EEG);
    ECG=cat(2,ecg_orig(1:fs(f_idx).*time_trim)',ECG);
    
    
    HeaderFile=headerFile{f_idx};
    MarkerFile=markerFile{f_idx};
    sfreq=fs(f_idx);

    EEG_orig=eeg_orig;
    ECG_orig=ecg_orig(:)';
        
    save(sprintf('%s.mat',fstem),'label','sfreq','HeaderFile','MarkerFile','EEG','ECG','EEG_orig','ECG_orig','time_trim','TRIGGER_ECG','TRIGGER','-v7.3','EEG_before_aas','EEG_before_bcg','EEG_bcg_pred','EEG_after_reref','EEG_after_bcg');
    

end;


return;

