close all; clear all;

headerFile = {
    '../eeg_raw/SSVEP_sms_1.vhdr';
    '../eeg_raw/SSVEP_sms_2.vhdr';
    '../eeg_raw/SSVEP_sms_3.vhdr';
    };

markerFile={
    '../eeg_raw/SSVEP_sms_1.vmrk';
    '../eeg_raw/SSVEP_sms_2.vmrk';
    '../eeg_raw/SSVEP_sms_3.vmrk';
    };
erp_event={1e3, [1 10]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG setup
%

TR=2; %second

%these two tokens are required (but values are arbitrary) for data
%collected inside MRI
trigger_token=1e3;
sync_token=1e2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f_idx=1:length(headerFile)
    [dummy,fstem]=fileparts(headerFile{f_idx});
    fprintf('reading [%s]...\n',fstem);


    % meta information such as samplingRate (fs), labels, etc
    [fs(f_idx) label meta] = bva_readheader(headerFile{f_idx});



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

    for event_idx=1:length(erp_event)
        tmp=erp_event{event_idx};
        trials=[];
        for ii=1:length(tmp)
            trials=union(trials,find(trigger{f_idx}.event==tmp(ii)));
        end;
        fprintf('\t[%d] events found for trigger [%s]...\n',length(trials),num2str(erp_event{event_idx}));
    end;


    file_para=sprintf('%s_soa.para',fstem);

    fprintf('writing [%s]....\n',file_para);
    fp=fopen(file_para,'w');

    soa_offset=0; %adjustment for slice-timing correction; second
    mri_idx=find(trigger{f_idx}.event==1e3);
    soa_start=min(trigger{f_idx}.time(mri_idx)./fs(f_idx));

    for idx=1:length(trigger{f_idx}.time)
        fprintf(fp,'%2.2f\t%d\n',trigger{f_idx}.time(idx)./fs(f_idx)+soa_offset-soa_start,trigger{f_idx}.event(idx));
        %fprintf(fp,'%2.2f\t%d\n',round((trigger{f_idx}.time(idx)./fs(f_idx)+soa_offset-soa_start)/TR)*TR,trigger{f_idx}.event(idx));
    end;

    fclose(fp);

end;


return;
