close all; clear all;

data_KIT={
    '../MEG_data/EXP_1.con';
    '../MEG_data/EXP_2.con';
    '../MEG_data/EXP_3.con';
    '../MEG_data/EXP_4.con';
    };

ft_defaults;

trigindx=[193 194];
meg_chan_idx=[1:157];

t_pre=0.1; %pre-stimulus interval (s)
t_post=1.0; %post-stimulus interval (s)

flag_erf_keep_raw_trials=0; %keep all trials

output_stem='erf_041019.mat';

for d_idx=1:length(data_KIT)
    
    fprintf('reading [%s]...\n',data_KIT{d_idx});
    hdr=ft_read_header(data_KIT{d_idx});
    dat = ft_read_data(data_KIT{d_idx}, 'header', hdr);
    [event] = ft_read_event(data_KIT{d_idx},'trigindx',trigindx);
    
    t_pre_samp=round(hdr.Fs*t_pre);
    t_post_samp=round(hdr.Fs*t_post);
    
    for trig_idx=1:length(trigindx)
        trig_str=sprintf('TRIG%d',trigindx(trig_idx));
        ev_list{trig_idx}=[];
        for ev_idx=1:length(event)
            if(strcmp(event(ev_idx).type, trig_str))
                ev_list{trig_idx}=cat(1,ev_list{trig_idx},ev_idx);
            end;
        end;
        
        
        fprintf('\t[%d] events found for trigger [%s]...\n',length(ev_list{trig_idx}),trig_str);
        
        
        fprintf('\tcollecting trials');
        for ev_idx=1:length(ev_list{trig_idx})
            if((event(ev_list{trig_idx}(ev_idx)).sample-t_pre_samp>0)&&(event(ev_list{trig_idx}(ev_idx)).sample+t_post_samp<size(dat,2)))
                fprintf('.');
                erf(d_idx,trig_idx).trial(:,:,ev_idx)=dat(:,event(ev_list{trig_idx}(ev_idx)).sample-t_pre_samp:event(ev_list{trig_idx}(ev_idx)).sample+t_post_samp);
            end;
        end;
        fprintf('\n');
        
        %calculate ERF
        erf(d_idx,trig_idx).erf=mean(erf(d_idx,trig_idx).trial(meg_chan_idx,:,:),3);
        
        %for ERF across data files
        if(d_idx==1) 
            erf_all(trig_idx).data=[];
        end;
        erf_all(trig_idx).data=cat(3,erf_all(trig_idx).data,erf(d_idx,trig_idx).trial(meg_chan_idx,:,:));
        
        %for baseline covariance matrix estimation
        if(d_idx==1) 
            erf_c_all(trig_idx).data=[];
        end;
        erf_c_all(trig_idx).data=cat(3,erf_c_all(trig_idx).data,erf(d_idx,trig_idx).trial(meg_chan_idx,1:t_pre_samp,:));
        
        %collect ERF info
        if(~flag_erf_keep_raw_trials)
            erf(d_idx,trig_idx).trial=[];
        end;
        
        erf(d_idx,trig_idx).timeVec=([1:size(erf(d_idx,trig_idx).erf,2)]-1)./hdr.Fs-t_pre;
        erf(d_idx,trig_idx).trig=trigindx(trig_idx);
        erf(d_idx,trig_idx).trig_str=trig_str;
        erf(d_idx,trig_idx).epoch_idx=ev_list{trig_idx};
        erf(d_idx,trig_idx).data_file=data_KIT{d_idx};
        
        erf_all(trig_idx).timeVec=([1:size(erf(d_idx,trig_idx).erf,2)]-1)./hdr.Fs-t_pre;
        erf_all(trig_idx).trig=trigindx(trig_idx);
        erf_all(trig_idx).trig_str=trig_str;
        erf_all(trig_idx).epoch_idx=ev_list{trig_idx};
        erf_all(trig_idx).data_file=data_KIT{d_idx};        
    end;
end;

for trig_idx=1:length(trigindx)
    erf_all(trig_idx).erf=mean(erf_all(trig_idx).data,3);
    erf_all(trig_idx).data=[];
    
    tmp=erf_c_all(trig_idx).data;
    tmp=reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)]);
    C(trig_idx).C=cov(tmp');
end;

save(output_stem,'erf','erf_all','C');