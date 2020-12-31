close all; clear all

%erf
file_erf='erf_050719.mat';

%fwd
file_fwd='seeg_fwd_wb_091019.mat';

output_stem='seeg_wb_mne_lambda_cv_123020';

SNR=[0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000];
%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);

%search correspondence
for d_idx=1:length(erf_all(1).name)
    n=erf_all(1).name{d_idx};
    n=n(2:end);
    n=erase(n,"'");
    
    %find the matched lead fields in the forward solution
    A_idx(d_idx)=find(strcmp(A(1).name, n));
end;

A2D=[];
for hemi_idx=1:2
    %make sure the lead field matched the data
    A(hemi_idx).A=A(hemi_idx).A(A_idx,:);
    
    %A(hemi_idx).A=A(hemi_idx).A(:,1:length(A(hemi_idx).v_idx)*3);
    
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

A2D_orig=A2D;

%noise covariance matrix
for idx=1:length(erf_all)
    if(idx==1) Ctmp=double(C(idx).C); else Ctmp=Ctmp+double(C(idx).C); end;
end;
C=Ctmp./length(erf_all);

C=diag(diag(C));
C_orig=C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrode information

file_mat='electrode_040119_085342.mat';
load(file_mat);
count=1;

for e_idx=1:length(electrode)
    
    for c_idx=1:electrode(e_idx).n_contact
        
        electrode_name=sprintf('%s%d',electrode(e_idx).name,c_idx);
        
        %if(~isempty(find(strcmp(electrode_name,ncov_electrode))))
        surface_coord(count,:)=electrode(e_idx).coord(c_idx,:)';
        count=count+1;
        %else
        %end;
    end;
end;
surface_coord=surface_coord(A_idx,:);
for c_idx=1:size(surface_coord,1)
    contact_dist_matrix(c_idx,:)=sqrt(sum((repmat(surface_coord(c_idx,:),[size(surface_coord,1),1])-surface_coord).^2,2));
end;

cv_select=[1:size(A2D_orig,1)]; %the index for the electrode to be left out

all_idx=[1:size(A2D_orig,1)];

data_select=[size(A2D_orig,1)-1, round((size(A2D_orig,1)).*0.9), round((size(A2D_orig,1)).*0.7), round((size(A2D_orig,1)).*0.5)]; %number of contacts used for source modeling

for snr_idx=1:length(SNR)
%     for data_idx=1:length(data_select)
%         
%         for cv_idx=1:length(cv_select)
%             
%             fprintf('SNR=%1.1e; Using [%d] contacts for source modeling; cross-validating contact [%03d]|[%03d]....\r', SNR(snr_idx),data_select(data_idx), cv_idx, length(cv_select));
%             
%             A2D=A2D_orig;
%             C=C_orig;
%             
%             %forward solution for the left-out contact
%             leftout_select=cv_select(cv_idx);
%             A2D_leftout=A2D(leftout_select,:);
%             
%             %forward solution for the contacts for source modeling
%             model_select=setdiff(all_idx,leftout_select); %all contacts other than the left-out contact
%             model_select=model_select(randperm(length(model_select))); %randomize the chosen contacts for source modeling
%             model_select=model_select(1:data_select(data_idx)); %select the first few contacts for source modeling
%            A2D=A2D_orig(model_select,:);
            A2D=A2D_orig;
            
%             contact_dist{data_idx}(cv_idx,:)=contact_dist_matrix(leftout_select,model_select);
%             contact_leftout_select{data_idx}(cv_idx,:)=leftout_select(:);
%             contact_model_select{data_idx}(cv_idx,:)=model_select(:);
            
            R2D=[];
            for hemi_idx=1:2
                R=ones(n_source(hemi_idx),1);
                %W=RA'*inv(A*R*A'+lambda.*C);
                RR=repmat(R,[1, 3])';
                R2D=cat(1,R2D,RR(:));
            end;
            
%            RAt=repmat(R2D,[1 length(model_select)]).*(A2D');
             RAt=repmat(R2D,[1 size(A2D,1)]).*(A2D');
            
            ARAt=A2D*RAt;
            
            p_signal=sum(diag(ARAt));
            
            %C=C(model_select,:);
            %C=C(:,model_select);
            p_noise=sum(diag(C));
            
            lambda=p_signal/p_noise/SNR(snr_idx);
            
            %gcv(snr_idx)=1./n.*norm(eye-A)./(
            
            W2D=RAt*inv(ARAt+lambda.*C);
            
%             %for trig_idx=1:size(erf_all,2)
             for trig_idx=1:1
                 Y=erf_all(trig_idx).erf;
                 baseline=find(erf_all(trig_idx).timeVec<0);
                 Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
%                 
%                 Y_leftout{snr_idx,data_idx}(cv_idx,:)=Y(leftout_select,:); %<-----left-out measurements
%                 
%                 
%                 X_mne0=W2D*Y(model_select,:);
                 A_lambda=eye(size(A2D,1))-A2D*W2D;
                 res=A_lambda*Y;
                 gcv_num=sum(res(:).^2)./size(Y,1);
                 gcv_den=(1./size(Y,1).*trace(A_lambda)).^2;
                 gcv(snr_idx)=gcv_num./gcv_den;
%                 Y_pred{snr_idx, data_idx}(cv_idx,:)=A2D_leftout*X_mne0; %<----predicted measurements
%                 
             end;
%             
%         end;
%     end;
end;


save wb_mne_lambda_gcv_123020.mat gcv SNR
return;

