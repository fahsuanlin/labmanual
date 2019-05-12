close all; clear all;

load erp_avg_outside.mat

tmp=[];
for i=1:size(erp,1)
     if(isempty(tmp))
        tmp=erp{i,3}.erp;
    else
        tmp=tmp+erp{i,3}.erp;
    end;
end;
data=tmp./size(erp,1);
timeVec=erp{1,3}.timeVec;

label=erp{1,3}.electrode_name;
flag_baseline_correct=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('eeg_fwd_prep_042319.mat','points');

%rotate the EEG channel position 90-deg at the x-y plane
T=[0 1 0 0 
  -1 0 0 0
   0 0 1 0
   0 0 0 1];
pt=[points ones(size(points,1),1)]';
pt=(T*pt).';
pt=pt(:,1:3);

hdr.grad.coilpos=pt.*1e3;

if(flag_baseline_correct)
    base_idx=find(timeVec<0);
    data=bsxfun(@minus,data,mean(data(:,base_idx),2));
end;

h=etc_plotEF_kit(data,'timeVec',timeVec,'label',label,'hdr',hdr,'threshold',[-5 5]);

return;



