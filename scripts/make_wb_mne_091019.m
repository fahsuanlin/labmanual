close all; clear all

%erf
file_erf='erf_050719.mat';
flag_auto_remove_comma=0; %remove ' from the electrode name automatically

%fwd
file_fwd='seeg_fwd_wb_091019.mat';

output_stem='seeg_wb_mne_091019';

SNR=10;
%%%%%%%%%%%%%%

load(file_erf);

load(file_fwd);

data_idx_remove=[];
%search correspondence
A_idx_tmp=zeros(1, length(erf_all(1).name)).*nan;
for d_idx=1:length(erf_all(1).name)
    n=erf_all(1).name{d_idx};
    n=n(2:end);
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
    Y=erf_all(trig_idx).erf;
    baseline=find(erf_all(trig_idx).timeVec<0);
    Y=bsxfun(@minus,Y,mean(Y(:,baseline),2));
    
    Y_pred=zeros(size(Y));

    Y(data_idx_remove,:)=[];

    X_mne0=W2D*Y;
   
    X_mne=reshape(X_mne0,[3 sum(n_source) size(Y,2)]);
    X_mne=squeeze(sqrt(sum(X_mne.^2,1)));
    clear X_mne0;
    
    X_dspm=bsxfun(@minus,X_mne, mean(X_mne(:,baseline),2));
    X_dspm=bsxfun(@rdivide,X_dspm,std(X_dspm(:,baseline),0,2));
    
    t0=mean(diff(erf_all(trig_idx).timeVec)); %sampling interval; ms;
    for hemi_idx=1:2
        switch hemi_idx
            case 1
                fn=sprintf('%s_%s-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(1:length(A(hemi_idx).v_idx),:);
                fn_mne=sprintf('%s_%s_mne-lh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi_mne=X_mne(1:length(A(hemi_idx).v_idx),:);
            case 2
                fn=sprintf('%s_%s-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi=X_dspm(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx),:);
                fn_mne=sprintf('%s_%s_mne-rh.stc',output_stem,erf_all(trig_idx).trig_str);
                X_hemi_mne=X_mne(n_source(1)+1:n_source(1)+length(A(hemi_idx).v_idx),:);
        end;
        fprintf('\tsaving [%s]...\n',fn);
        inverse_write_stc(X_hemi,A(hemi_idx).v_idx,min(erf_all(trig_idx).timeVec),t0,fn);
        fprintf('\tsaving [%s]...\n',fn_mne);
        inverse_write_stc(X_hemi_mne,A(hemi_idx).v_idx,min(erf_all(trig_idx).timeVec),t0,fn_mne);
    end;
    
    fn=sprintf('%s_%s-vol.stc',output_stem,erf_all(trig_idx).trig_str);
    fprintf('\tsaving [%s]...\n',fn);
    inverse_write_stc(X_dspm,[0:size(X_dspm,1)-1],min(erf_all(trig_idx).timeVec),t0,fn);
    
    
%     subjects_dir=getenv('SUBJECTS_DIR');
%     subject='s036';
%     surf='orig';
%     overlay_smooth=5;
%     
%     file_source_wholebrain='/Users/fhlin_admin/workspace/seeg/subjects/s036/mri/orig.mgz';
%     src_wb=MRIread(file_source_wholebrain);
% 
% 
%     %initialize
%     loc_vol=[];
%     for hemi_idx=1:2
%         switch hemi_idx
%             case 1
%                 offset=0;
%                 hemi='lh';
%             case 2
%                 offset=n_source(1);
%                 hemi='rh';
%         end;
%         
%         for time_idx=1:size(X_dspm,2)
%             %smooth overaly over cortical surface
%             if(time_idx==1)
%                 surf='orig';
%                 file_surf=sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,hemi,surf);
%                 [vertex_coords, faces] = read_surf(file_surf);
%                 
%                 overlay_D=[];
%             end;
%             
%             if(mod(time_idx,50)==0) fprintf('[%1.1f%%]...\r',time_idx./size(X_dspm,2).*100); end;
%             
%             X_hemi_cort=X_dspm(offset+1:offset+length(A(hemi_idx).v_idx),time_idx);
%             X_hemi_subcort=X_dspm(offset+length(A(hemi_idx).v_idx)+1:offset+n_source(hemi_idx),time_idx);
%             ov=zeros(size(vertex_coords,1),1);
%             ov(A(hemi_idx).v_idx+1)=X_hemi_cort;
%             
%             [ovs,dd0,dd1,overlay_Ds,overlay_D]=inverse_smooth('','vertex',vertex_coords','face',faces','value',ov,'value_idx',A(hemi_idx).v_idx+1,'step',overlay_smooth,'n_ratio',length(ov)/size(X_hemi_cort,1),'flag_display',0,'flag_regrid',0,'flag_fixval',0,'D',overlay_D);
%             
%             %collect "smoothed" overlay values at cortical and subortical areas
%             %X_wb=cat(1,X_wb,cat(1,ovs(:),X_hemi_subcort(:)));
%             if(time_idx==1)
%                 X_wb{hemi_idx}=zeros(length(ovs(:))+length(X_hemi_subcort(:)),size(X_dspm,2));
%             end;
%             X_wb{hemi_idx}(:,time_idx)=cat(1,ovs(:),X_hemi_subcort(:));
%             
%             if(time_idx==1)
%                 %get coordinates from surface to volume
%                 loc=cat(1,vertex_coords./1e3,A(hemi_idx).wb_loc);
%                 loc_surf=[loc.*1e3 ones(size(loc,1),1)]';
%                 tmp=inv(src_wb.tkrvox2ras)*loc_surf;
%                 %loc_vol=cat(1,loc_vol,round(tmp(1:3,:))');
%                 loc_vol{hemi_idx}=round(tmp(1:3,:))';
%                 
%                 loc_vol_idx{hemi_idx}=sub2ind(size(src_wb.vol),loc_vol{hemi_idx}(:,2),loc_vol{hemi_idx}(:,1),loc_vol{hemi_idx}(:,3));
% 
%             end;
%         end;
%         fprintf('\n');
%     end;
%     
%     for time_idx=1:size(X_dspm,2)
%         if(mod(time_idx,50)==0) fprintf('[%1.1f%%]***\r',time_idx./size(X_dspm,2).*100); end;
%         
%         %show overlay with T1 underlay
%         if(time_idx==1)
%             sz=size(src_wb.vol);
%             overlay_vol=zeros(sz(1),sz(2),sz(3),size(X_dspm,2));
%         end;
%         tmp=zeros(size(src_wb.vol));
%         
%         for hemi_idx=1:2
%             tmp(loc_vol_idx{hemi_idx})=X_wb{hemi_idx}(:,time_idx);
%         end;
%         
%         overlay_vol(:,:,:,time_idx)=tmp;
%         overlay_vol(:,:,:,time_idx)=imgaussfilt3(overlay_vol(:,:,:,time_idx),2);
%     end;
%     fprintf('\n');
    
end;
