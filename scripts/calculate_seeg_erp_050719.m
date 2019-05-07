close all; clear all;

da=dir('../epochs/B1(10)_epoch*.mat'); %audio
dv=dir('../epochs/B2(20)_epoch*.mat'); %visual
dav=dir('../epochs/B3(30)_epoch*.mat'); %audiovisual

%contacts info: ../epochs/seeg_channels.txt
[dummy, contact_name]=textread('../epochs/seeg_channels.txt','%d %s','headerlines',1);

%time points info: ../epochs/seeg_epoch_timepoints.txt
[dummy, timeVec]=textread('../epochs/seeg_epoch_timepoints.txt','%d %f','headerlines',1);

T_contact_idx=[35:40]; %sT1' to sT6'
mont_config=[2 4];

file_output='erf_050719';

flag_reref=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx=1:length(da)
    fprintf('loading [%s]...\r',da(idx).name);
    load(sprintf('../epochs/%s',da(idx).name));
    if(flag_reref)
        %re-reference
        epoch=bsxfun(@minus,epoch,mean(epoch,1));
    end;
    epoch_a(:,:,idx)=epoch;
end;
fprintf('\n');
epoch_a_avg=mean(epoch_a,3);
epoch_a_std=std(epoch_a,0,3);

for idx=1:length(dv)
    fprintf('loading [%s]...\r',dv(idx).name);
    load(sprintf('../epochs/%s',dv(idx).name));
    if(flag_reref)
        %re-reference
        epoch=bsxfun(@minus,epoch,mean(epoch,1));
    end;
    epoch_v(:,:,idx)=epoch;
end;
fprintf('\n');
epoch_v_avg=mean(epoch_v,3);
epoch_v_std=std(epoch_v,0,3);

for idx=1:length(dav)
    fprintf('loading [%s]...\r',dav(idx).name);
    load(sprintf('../epochs/%s',dav(idx).name));
    if(flag_reref)
        %re-reference
        epoch=bsxfun(@minus,epoch,mean(epoch,1));
    end;
    epoch_av(:,:,idx)=epoch;
end;
fprintf('\n');
epoch_av_avg=mean(epoch_av,3);
epoch_av_std=std(epoch_av,0,3);


epoch_av_a_avg=epoch_av_avg-epoch_a_avg;
epoch_av_v_avg=epoch_av_avg-epoch_v_avg;


erf_all(1).erf=epoch_a_avg;
erf_all(2).erf=epoch_v_avg;
erf_all(3).erf=epoch_av_avg;

erf_all(1).name=contact_name;
erf_all(2).name=contact_name;
erf_all(3).name=contact_name;

erf_all(1).trig_str='a';
erf_all(2).trig_str='v';
erf_all(3).trig_str='av';


erf_all(1).timeVec=timeVec;
erf_all(2).timeVec=timeVec;
erf_all(3).timeVec=timeVec;

baseline_idx=find(erf_all(1).timeVec<0);
tmp=epoch_a(:,baseline_idx,:);
tmp=reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)])';
C(1).C=cov(tmp);
tmp=epoch_v(:,baseline_idx,:);
tmp=reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)])';
C(2).C=cov(tmp);
tmp=epoch_av(:,baseline_idx,:);
tmp=reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)])';
C(3).C=cov(tmp);

save(file_output,'erf_all','C');
return;
