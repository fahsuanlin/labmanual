close all; clear all;

data_file_filter={
	'../epoch_mat/sH''1-sH''2-1.0*.mat';
        '../epoch_mat/sH''1-sH''2-2.0*.mat';
        '../epoch_mat/sH''1-sH''2-3.0*.mat';
};

data_name={
	'sH1-sH2-1p0';
        'sH1-sH2-2p0';
        'sH1-sH2-3p0';
};


%contacts info: ../epochs/seeg_channels.txt
[dummy, contact_name]=textread('../epoch_mat/seeg_channels.txt','%d %s','headerlines',1);

%time points info: ../epochs/seeg_epoch_timepoints.txt
[dummy, timeVec]=textread('../epoch_mat/seeg_epoch_timepoints.txt','%d %f','headerlines',1);

%T_contact_idx=[35:40]; %sT1' to sT6'
%mont_config=[2 4];

file_output='erf_raw_042121';

flag_reref=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for filter_idx=1:length(data_file_filter)
	dd=dir(data_file_filter{filter_idx});
	for idx=1:length(dd)
    		fprintf('loading [%s]...\r',dd(idx).name);
    		load(sprintf('%s/%s',fileparts(data_file_filter{filter_idx}),dd(idx).name));
    		if(flag_reref)
        		%re-reference
        		epoch=bsxfun(@minus,epoch,mean(epoch,1));
    		end;
    		epoch(:,:,idx)=epoch;
	end;
	fprintf('\n');
	epoch_avg=mean(epoch,3);
	epoch_std=std(epoch,0,3);

	erf_all(filter_idx).erf=epoch_avg;
	erf_all(filter_idx).name=contact_name;
	erf_all(filter_idx).trig_str=data_name{filter_idx};
	erf_all(filter_idx).erf_raw=epoch;
	erf_all(filter_idx).timeVec=timeVec;

	baseline_idx=find(erf_all(filter_idx).timeVec<0);
	tmp=epoch(:,baseline_idx,:);
	tmp=reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)])';
	C(filter_idx).C=cov(tmp);
	C(filter_idx).name=contact_name;
end;

save(file_output,'erf_all','C');
return;

T_idx=find(not(cellfun('isempty',strfind(erf_all(1).name,'T'))));
for e_idx=1:length(erf_all)
	erf_T(e_idx)=erf_all(e_idx);
	erf_T(e_idx).name=contact_name(T_idx);
	erf_T(e_idx).erf=erf_all(e_idx).erf(T_idx,:);
	erf_T(e_idx).erf_raw=erf_all(e_idx).erf_raw(T_idx,:,:);

	C_T(e_idx)=C(e_idx);
	C_T(e_idx).C=C_T(e_idx).C(T_idx,T_idx);
	C_T(e_idx).name=C_T(e_idx).name(T_idx);
end;


save(sprintf('%s_T',file_output),'erf_T','C_T');

return;
