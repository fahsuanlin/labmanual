close all; clear all

%erf
file_erf='erf_061721.mat';
flag_auto_remove_comma=0; %remove ' from the electrode name automatically

%fwd
file_fwd='seeg_fwd_wb_dec_061721.mat';

output_stem='seeg_wb_mne_cc_roi_raw_061721';

SNR=100;


file_label={
    'aud_HCPMMP1-lh_s031.label';
    'aud_a1_HCPMMP1-lh_s031.label';
    'aud_belt_HCPMMP1-lh_s031.label';
    'aud_HCPMMP1-rh_s031.label';
    'aud_a1_HCPMMP1-rh_s031.label';
    'aud_belt_HCPMMP1-rh_s031.label';
    'lh.G_temp_sup-G_T_transv_s031.label';
    'lh.G_temp_sup-Lateral_div1_s031.label';
    'lh.G_temp_sup-Lateral_div2_s031.label';
    'lh.planum_temporale_s031.label';
    'lh.S_intrapariet_and_P_trans_s031.label';
    'lh.S_occipital_ant_s031.label';
    'lh.S_temporal_sup_div1_s031.label';
    'lh.S_temporal_sup_div2_s031.label';
    'lh.S_temporal_sup_div3_s031.label';
    'lh.S_temporal_transverse_s031.label';
    'rh.G_temp_sup-G_T_transv_s031.label';
    'rh.G_temp_sup-Lateral_div1_s031.label';
    'rh.G_temp_sup-Lateral_div2_s031.label';
    'rh.planum_temporale_s031.label';
    'rh.S_intrapariet_and_P_trans_s031.label';
    'rh.S_occipital_ant_s031.label';
    'rh.S_temporal_sup_div1_s031.label';
    'rh.S_temporal_sup_div2_s031.label';
    'rh.S_temporal_sup_div3_s031.label';
    'rh.S_temporal_transverse_s031.label';
};


%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);

data_idx_remove=[];


for label_idx=1:length(file_label)
	[dummy,fstem]=fileparts(file_label{label_idx});
	fprintf('loadin <%s>...',file_label{label_idx});
	[roi(label_idx).idx]=inverse_read_label(file_label{label_idx},'flag_display',0);
	roi(label_idx).name=fstem;
	if(isempty(findstr(fstem,'lh.'))&&isempty(findstr(fstem,'-lh')))
		fprintf('[RH]...\n');
		roi(label_idx).hemi='rh';
	else
                fprintf('[LH]...\n');
                roi(label_idx).hemi='lh';
	end;
end;


%search correspondence
A_idx_tmp=zeros(1, length(erf_all(1).name)).*nan;
for d_idx=1:length(erf_all(1).name)
    n=erf_all(1).name{d_idx};
    %n=n(2:end);
    if(flag_auto_remove_comma)
	    n=erase(n,a);
    end;

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
    %W=RA'*inv(A*R*A'+lambda.*C);
    RR=repmat(R,[1, 3])';
    R2D=cat(1,R2D,RR(:));
end;

RAt=repmat(R2D,[1 n_chan]).*(A2D');

ARAt=A2D*RAt;

p_signal=sum(diag(ARAt));

%noise covariance matrix
for idx=1:length(erf_all)
    if(idx==1) Ctmp=C(idx).C; else Ctmp=Ctmp+C(idx).C; end;
end;
C=Ctmp./length(erf_all);
C(data_idx_remove,:)=[];
C(:,data_idx_remove)=[];

C=diag(diag(C));
p_noise=sum(diag(C));

lambda=p_signal/p_noise/SNR;

W2D=RAt*inv(ARAt+lambda.*C);

for roi_idx=1:length(roi)
	if(strcmp(roi(roi_idx).hemi,'lh'))
		[dummy,vv]=intersect(A(1).v_idx,roi(roi_idx).idx);
		offset=0;
                vv1=offset+(vv-1).*3+1;
                vv2=offset+(vv-1).*3+2;
                vv3=offset+(vv-1).*3+3;
                vv_cat=cat(2,vv1,vv2,vv3)';
                roi(roi_idx).W=W2D(vv_cat(:),:);
	else
		[dummy,vv]=intersect(A(2).v_idx,roi(roi_idx).idx);
		offset=size(A(1).A,2);
                vv1=offset+(vv-1).*3+1;
                vv2=offset+(vv-1).*3+2;
                vv3=offset+(vv-1).*3+3;
                vv_cat=cat(2,vv1,vv2,vv3)';
	        roi(roi_idx).W=W2D(vv_cat(:),:);
	end;
end;




for trig_idx=1:size(erf_all,2)
    for trial_idx=1:size(erf_all(trig_idx).erf_raw,3)
	fprintf('trigger [%s]:: <%04d|%04d>...\r',erf_all(trig_idx).trig_str,trial_idx,size(erf_all(trig_idx).erf_raw,3));
	Y=erf_all(trig_idx).erf_raw(:,:,trial_idx);
    	%baseline=find(erf_all(trig_idx).timeVec<0);
    	%Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));

    	Y(data_idx_remove,:)=[];

    	for roi_idx=1:length(roi)
                mne=roi(roi_idx).W*Y;
                mne=reshape(mne,[3,size(mne,1)./3,size(mne,2)]);
                mne=squeeze(mean(mne,2));

                roi(roi_idx).X_mne{trig_idx}(:,trial_idx,:)=mne;

                roi(roi_idx).X_mne_rms{trig_idx}(trial_idx,:)=squeeze(sqrt(sum(roi(roi_idx).X_mne{trig_idx}(:,trial_idx,:).^2,1)));
                roi(roi_idx).timeVec{trig_idx}=erf_all(trig_idx).timeVec;
                roi(roi_idx).condition{trig_idx}=erf_all(trig_idx).trig_str;
    	end;
    end;
	fprintf('\n');
end;

save(output_stem,'-v7.3','roi');

