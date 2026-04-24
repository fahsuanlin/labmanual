close all; clear all

%erf
file_erf='erf_061721.mat';
flag_auto_remove_comma=0; %remove ' from the electrode name automatically

%fwd
file_fwd='seeg_fwd_wb_dec_061721.mat';

output_stem='seeg_wb_mne_cc_raw_061721';

subject='s031';
target_subject='fsaverage';

SNR=100;
%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);

data_idx_remove=[];
%search correspondence
A_idx_tmp=zeros(1, length(erf_all(1).name)).*nan;
for d_idx=1:length(erf_all(1).name)
    n=erf_all(1).name{d_idx};
    %n=n(2:end);
    if(flag_auto_remove_comma)
	    n=erase(n,"'");
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

for trig_idx=1:size(erf_all,2)
	for dd=1:3 %x,y,z
		switch dd
		case 1
			dd_str='x';
		case 2
			dd_str='y';
		case 3
			dd_str='z';
		end;
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    fn_mne=sprintf('%s_%s_mne_%s-lh.stc',output_stem,erf_all(trig_idx).trig_str,dd_str);
                    hemi='lh';
                case 2
                    fn_mne=sprintf('%s_%s_mne_%s-rh.stc',output_stem,erf_all(trig_idx).trig_str,dd_str);
                    hemi='rh';
            end;
            %fprintf('\tsaving [%s]...\n',fn);
            fprintf('\tloading [%s]...\n',fn_mne);
            [stc{hemi_idx}(:,:,dd),v{hemi_idx},a,b,c]=inverse_read_stc(fn_mne);
        end;
        end;
        

	for hemi_idx=1:2
		for t_idx=1:size(stc{hemi_idx},2)
			stc_cort{hemi_idx}(:,t_idx)=sum(squeeze(stc{hemi_idx}(:,t_idx,:)).*A(hemi_idx).ori,2);
		end;
		switch hemi_idx
                case 1
                    fn_mne=sprintf('%s_%s_mne_cort-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                    hemi='lh';
                case 2
                    fn_mne=sprintf('%s_%s_mne_cort-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                    hemi='rh';
            	end;
		fprintf('\tsaving [%s]...\n',fn_mne);
		inverse_write_stc(stc_cort{hemi_idx},v{hemi_idx},a,b,fn_mne);
	end;

        %morphing STC....
		for hemi_idx=1:2
            		switch hemi_idx
                	case 1
                   		hemi='lh';
                	case 2
                	   	hemi='rh';
            		end;
			fn_in=sprintf('%s_%s_mne_cort-%s.stc',output_stem,erf_all(trig_idx).trig_str,hemi);
            		fn_out=sprintf('%s_2_%s_%s_%s_mne_cort',subject,target_subject,output_stem,erf_all(trig_idx).trig_str);
            		cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject, fn_in, target_subject, fn_out, hemi);
            		eval(cmd);
        	end;
%    end;
    
end;
