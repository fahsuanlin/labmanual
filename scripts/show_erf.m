close all; clear all;

data_KIT={
    '../MEG_data/EXP_1.con';
    };
hdr=ft_read_header(data_KIT{1});

load('erf_041019.mat');

data=erf_all(1).erf.*1e15; %fT
timeVec=erf_all(1).timeVec*1e3; %ms
label=hdr.grad.label;
flag_baseline_correct=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(flag_baseline_correct)
    base_idx=find(timeVec<0);
    data=bsxfun(@minus,data,mean(data(:,base_idx),2));
end;

h=etc_plotEF_kit(data,'timeVec',timeVec,'label',label);

return;



