close all; clear all;


file_eeg={
    'SSVEP_noMR_1_ccm_042124.mat';
    };

file_output='erp_inside_ccm.mat';
file_avg_output='erp_avg_inside_ccm.mat';

select_channel={ 'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'    'ECG'};

erp_pre=0.2; %s; pre-stimulus interval
erp_post=1.0; %s; post-stimulus interval
flag_badrejection=1; %automatic bad trial rejection
badrejection_threshold=100; %microV; threshold to consider as a bad trial
flag_baseline_corr=0;
erp_event={2, 3,[2 3]};

for f_idx=1:length(file_eeg)
    load(file_eeg{f_idx});
    
    
    trigger{f_idx}=TRIGGER;
    fs(f_idx)=sfreq;
    eeg{f_idx}=EEG;
    
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
                    tmp=ones(size(tmp)).*inf;
                end;
                
                epoch_abs_max(:,epoch_idx)=max(abs(tmp),[],2);
            end
            reject_trial=find(max(epoch_abs_max,[],1)>badrejection_threshold);
            fprintf('trials {%s} with maximum higher than %1.1f (uV).\n',mat2str(reject_trial),badrejection_threshold);
            %epoch_data(:,:,reject_trial)=[];
            %trigger{f_idx}.event(reject_trial)=[];
        else
            reject_trial=[];
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
        if(~iscell(trigger{f_idx}.event))
            str={};
            for idx=1:length(trigger{f_idx}.event)
                str{idx}=sprintf('%d',trigger{f_idx}.event(idx));
            end;
            trigger{f_idx}.event=str;
        end;
        
        for event_idx=1:length(erp_event)
            tmp=erp_event{event_idx};
            str={}; for i=1:length(tmp) str{i}=sprintf('%d',tmp(i)); end; tmp=str;
            trials=[];
            for ii=1:length(tmp)
                %trials=union(trials,find(trigger{f_idx}.event==tmp(ii)));
                trials=union(trials,find(strcmp(trigger{f_idx}.event, tmp{ii})));
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

%save(file_output,'erp');
save(file_avg_output,'erp_avg');
