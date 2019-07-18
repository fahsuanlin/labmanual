close all; clear all;

headerFile = {
    '../eeg_raw/sms01.vhdr';
    '../eeg_raw/sms02.vhdr';
    '../eeg_raw/sms03.vhdr';
    '../eeg_raw/sms04.vhdr';
    };

markerFile={
    '../eeg_raw/sms01.vmrk';
    '../eeg_raw/sms02.vmrk';
    '../eeg_raw/sms03.vmrk';
    '../eeg_raw/sms04.vmrk';
    };

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

flag_reref= [1 1 1 1]; %referencing: remove the average time course across all electrodes from the time course at each individual electrode
flag_hp=    [1 1 1 1]; %high pass filtering at 0.1 Hz
flag_aas=   [1 1 1 1]; %AAS for gradient artifact suppression
flag_bcg=    [1 1 1 1]; %BCG artifact suppression

time_trim=[20];     %trim off the first few seconds...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERP setup
%
erp_pre=0.2; %s; pre-stimulus interval
erp_post=1.5; %s; post-stimulus interval
flag_badrejection=1; %automatic bad trial rejection
flag_baseline_corr=1;
erp_event={1,10, [1 10]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file archiving setup
%

%output_file='erp.mat';
output_avg_file='erp_avg_smsini.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f_idx=1:length(headerFile)
    [dummy,fstem]=fileparts(headerFile{f_idx});
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
    
    ecg{f_idx}=eeg{f_idx}(ecg_channel,:)';
    eeg{f_idx}=eeg{f_idx}(eeg_channel,:);
    
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
    
    %high-pass filtering
    if(flag_hp(f_idx))
        fprintf('\tHP...\n');
        %high-pass filtering (0.1 Hz)
        Wn = 0.1*2/fs(f_idx);
        N = 3; % order of 3 less processing
        [a,b] = butter(N,Wn,'high'); %bandpass filtering
        for ch_idx=1:size(eeg{f_idx},1)
            eeg{f_idx}(ch_idx,:) = filtfilt(a,b,eeg{f_idx}(ch_idx,:));
        end;
        if(~isempty(eeg_ref{f_idx}))
            eeg_ref{f_idx}=filtfilt(a,b,eeg_ref{f_idx}(:))';
        end;
        
        ecg{f_idx}=filtfilt(a,b,ecg{f_idx}(:))';
    end;
    
    
    %read maker file
    fprintf('\treading triggers...\n');
    mk=textread(markerFile{f_idx},'%s','headerlines',12,'delimiter','\n');
    if(~isempty(mk))
        for mk_idx=1:length(mk)
            tmp=strfind(mk{mk_idx},',');
            for tmp_idx=1:length(tmp)+1
                if(tmp_idx==1)
                    bb=1;
                else
                    bb=tmp(tmp_idx-1)+1;
                end;
                if(tmp_idx==length(tmp)+1)
                    ee=length(mk{mk_idx});
                else
                    ee=tmp(tmp_idx)-1;
                end;
                mk_tmp{tmp_idx}=mk{mk_idx}(bb:ee);
            end;
            
            if(strcmp(mk_tmp{2},'R128')) %MRI trigger
                trigger{f_idx}.event(mk_idx)=trigger_token;
                trigger{f_idx}.time(mk_idx)=str2num(mk_tmp{3}); %data index for the trigger
            elseif(~isempty(findstr(mk_tmp{2},'Sync'))) %sync
                trigger{f_idx}.event(mk_idx)=sync_token;
                trigger{f_idx}.time(mk_idx)=str2num(mk_tmp{3}); %data index for the trigger
            else %true events
                trigger{f_idx}.event(mk_idx)=str2num(mk_tmp{2}(2:end)); %trigger
                trigger{f_idx}.time(mk_idx)=str2num(mk_tmp{3}); %data index for the trigger
            end;
        end;
        dummy=sort(trigger{f_idx}.event);
        events=[dummy(find(diff(dummy))),dummy(end)];
        fprintf('\t\ttotal [%d] events found : {%s}\n',length(events),mat2str(events));
    else
        trigger{f_idx}=[];
    end;
    
    
    ch20_before_aas=eeg{f_idx}(20,:);
    ecg_before=ecg{f_idx};
    if(flag_aas(f_idx))
        fprintf('\tAAS...\n');
        gradient_trigger{f_idx}=zeros(size(eeg{f_idx},2),1);
        iidx=find(trigger{f_idx}.event==trigger_token);
        gradient_trigger{f_idx}(trigger{f_idx}.time(iidx))=1e3;
        ecg{f_idx}=eeg_ga(ecg{f_idx},gradient_trigger{f_idx},TR,fs(f_idx),'flag_display',0,'flag_ma_aas',1,'flag_aas_svd',0,'flag_anchor_bnd',0); %AAS on ECG
        ecg{f_idx}=sgolayfilt(ecg{f_idx},6,301); %smoothing out residual GA in ECG
        eeg{f_idx}=eeg_ga(eeg{f_idx},gradient_trigger{f_idx},TR,fs(f_idx),'flag_display',0,'flag_ma_aas',1,'flag_aas_svd',0,'flag_anchor_bnd',0); %AAS on EEG
    end;
    ch20_after_aas=eeg{f_idx}(20,:);
    ecg_after_aas=ecg{f_idx};
    
    %remove first few seconds (if needed....)
    if(~isempty(time_trim))
        time_trim_idx=round(time_trim*fs(f_idx));
        
        fprintf('trimming [%1.1f] s data {(%d) samples}....\n',time_trim, time_trim_idx);
        ecg{f_idx}=ecg{f_idx}(time_trim_idx+1:end);
        eeg{f_idx}=eeg{f_idx}(:,time_trim_idx+1:end);
        trigger{f_idx}.time=trigger{f_idx}.time-time_trim_idx;
        
        %update trigger info
        idx=find(trigger{f_idx}.time<0);
        trigger{f_idx}.time(idx)=[];
        trigger{f_idx}.event(idx)=[];
    end;

    [qrs_amp_raw,qrs{f_idx}]=pan_tompkin(ecg{f_idx},fs(f_idx),0);
    %%plot ECG and QRS detection
    %tt=[1:length(ecg{f_idx})]./fs(f_idx);
    %figure; plot(tt,ecg{f_idx}); hold on;
    %line(repmat(tt(qrs{f_idx}),[2 1]),repmat([min(ecg{f_idx}-mean(ecg{f_idx}))/2; max(ecg{f_idx}-mean(ecg{f_idx}))/2],size(qrs{f_idx})),'LineWidth',2.5,'LineStyle','-.','Color','r');
    
    if(flag_bcg(f_idx))
        fprintf('\tBCG...\n');
        %[eeg{f_idx}, qrs_i_raw]=eeg_bcg(eeg{f_idx},ecg{f_idx},fs(f_idx),'flag_display',0); %BCG correction on EEG
        [eeg{f_idx}, qrs_i_raw, bcg_all, ecg_all, bad_trials]=eeg_bcg(eeg{f_idx},ecg{f_idx},fs(f_idx),'flag_display',0,'bcg_nsvd',3,'n_ma_bcg',20,'flag_anchor_ends',0,'flag_post_ssp',0,'flag_badrejection',0); %BCG correction on EEG
    end;
    ch20_after_aas_bcg=eeg{f_idx}(20,:);
    
    
    trigger_ecg{f_idx}=trigger{f_idx};
    trigger_ecg{f_idx}.event=cat(2,trigger_ecg{f_idx}.event,33.*ones(1,length(qrs{f_idx}(:))));
    trigger_ecg{f_idx}.time=cat(2,trigger_ecg{f_idx}.time,qrs{f_idx}(:)');

    EEG=eeg{end};
    ECG=ecg{end};
    TRIGGER_ECG=trigger_ecg{end};
    TRIGGER=trigger{end};
    save(sprintf('%s.mat',fstem),'EEG','ECG','TRIGGER_ECG','TRIGGER','-v7.3');

    %epoching
    if(~isempty(trigger{f_idx}))
        time_pre=round(erp_pre.*fs(f_idx));
        time_total=round((erp_post+erp_pre).*fs(f_idx));
        epoch_data=zeros(size(eeg{f_idx},1),time_total,length(trigger{f_idx}.event)); %erp: channel x time x trials
        epoch_timeVec=([1:time_total]-1)./fs(f_idx)-erp_pre;
        fprintf('\tERP {%d} epoching',f_idx);
        for epoch_idx=1:length(trigger{f_idx}.event)
            fprintf('.');
            if((trigger{f_idx}.time(epoch_idx)-time_pre)>1)
                if(trigger{f_idx}.time(epoch_idx)+time_total-time_pre<size(eeg{f_idx},2))
                    epoch_data(:,:,epoch_idx)=eeg{f_idx}(:,trigger{f_idx}.time(epoch_idx)-time_pre:trigger{f_idx}.time(epoch_idx)+time_total-time_pre-1);
                 end;
            end;
        end;
        fprintf('\n');
        
        %automatic bad trial rejection
        if(flag_badrejection)
            fprintf('\tbad trial rejection...',f_idx);
            n_badtrial=0;
            for epoch_idx=1:length(trigger{f_idx}.event)
                tmp=epoch_data(eeg_channel,:,epoch_idx);
                
                if(~isempty(find(isnan(tmp))))
                    n_badtrial=n_badtrial+1;
                    tmp=ones(size(tmp)).*1e3;
                end;
                
                epoch_abs_max(:,epoch_idx)=max(abs(tmp),[],2);
            end
            reject_trial=find(max(epoch_abs_max,[],1)>700);
            fprintf('trials {%s} with maximum higher than 700 (uV).\n',mat2str(reject_trial));
            %epoch_data(:,:,reject_trial)=[];
            %trigger{f_idx}.event(reject_trial)=[];
        end;
        
        %baseline correction
        if(flag_baseline_corr)
            fprintf('\tbaseline correction',f_idx);
            baseline_idx=find(epoch_timeVec<0);
            for epoch_idx=1:size(epoch_data,3)
                fprintf('.');
                epoch_data(:,:,epoch_idx)=epoch_data(:,:,epoch_idx)-repmat(squeeze(mean(epoch_data(:,baseline_idx,epoch_idx),2)),[1,size(epoch_data,2)]);
            end
            fprintf('\n');
        end;
        
        %getting ERP
        for event_idx=1:length(erp_event)
            tmp=erp_event{event_idx};
            trials=[];
            for ii=1:length(tmp)
                trials=union(trials,find(trigger{f_idx}.event==tmp(ii)));
            end;
            trials=setdiff(trials,reject_trial);
            fprintf('\t[%d] events found for trigger [%s]...\n',length(trials),num2str(erp_event{event_idx}));
            erp{f_idx,event_idx}.trials=epoch_data(:,:,trials);
            erp{f_idx,event_idx}.n_trial=length(trials);
            erp{f_idx,event_idx}.trial_idx=trials;
            erp{f_idx,event_idx}.erp=mean(epoch_data(:,:,trials),3);
            erp{f_idx,event_idx}.timeVec=epoch_timeVec;
            erp{f_idx,event_idx}.trigger=erp_event(event_idx);
            erp{f_idx,event_idx}.electrode_name=found_channel;
            
            erp_avg{f_idx,event_idx}=erp{f_idx,event_idx};
            erp_avg{f_idx,event_idx}.trials=[];
        end;
    else
        erp=[];
    end;
end;

%save(output_file,'erp','fs','found_channel');

%save(output_avg_file,'erp','fs','found_channel');

return;

