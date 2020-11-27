close all; clear all;

fstem={
'isc_permswb_group_121419_z_hypothesis01_gmr_median_1mp';
'isc_permswb_group_121419_z_hypothesis02_gmr_median_1mp';
'isc_permswb_group_121419_z_hypothesis03_gmr_median_1mp';
'isc_permswb_group_121419_z_hypothesis04_gmr_median_1mp';
    };

threshold={
[0.95 0.99];
[0.95 0.99];
[0.95 0.99];
[0.95 0.99];
};

for f_idx=1:length(fstem)
    
    [stc_lh,v_lh]=inverse_read_stc(sprintf('%s-lh.stc',fstem{f_idx}));
    [stc_rh,v_rh]=inverse_read_stc(sprintf('%s-rh.stc',fstem{f_idx}));
    
    
    %etc_render_fsbrain('overlay_value',stc_lh(:,1),'overlay_vertex',v_lh,'overlay_threshold',[5 10]);
    %return;
    
    stc=cat(1,stc_lh(:,1),stc_rh(:,1));
    ii=find(isnan(stc(:)));
    stc(ii)=randn(size(ii)).*eps;
    
    idx=find(stc>0);
    p_uncorr(idx)=2-2*normcdf(stc(idx));
    idx=find(stc<0);
    p_uncorr(idx)=2-2*normcdf(-stc(idx));
    idx=find(p_uncorr<=0.05);
    
%     [pID_l]=FDR(p_uncorr(:),max(threshold{f_idx})./2);
%     Zid_l=norminv(1-pID_l/2);
%     v_l=length(find(abs(stc(:))>Zid_l));
%     
%     [pID_h]=FDR(p_uncorr(:),min(threshold{f_idx})./2);
%     Zid_h=norminv(1-pID_h/2);
%     v_h=length(find(abs(stc(:))>Zid_h));
    
%    fprintf('[%s]: [%d] voxel over threshold Z=%1.2f  (corr. p=%2.2f)\n',fstem{f_idx},v_l,Zid_l,max(threshold{f_idx}));
%    fprintf('[%s]: [%d] voxel over threshold Z=%1.2f  (corr. p=%2.2f)\n',fstem{f_idx},v_h,Zid_h,min(threshold{f_idx}));
    
    %     Zid_l=4.0;
    %     Zid_h=4.3;
    %
    %
        if(1)
            figure;
%           etc_render_fsbrain_stc({fstem{f_idx}},[Zid_l Zid_h],'overlay_smooth',5,'overlay_exclude_fstem','exclude','surf','inflated','flag_curv',1);
            etc_render_fsbrain_stc({fstem{f_idx}},threshold{f_idx},'overlay_smooth',5,'overlay_exclude_fstem','exclude','surf','inflated','flag_curv',1);
            set(gcf,'pos',[680   770   960   400]);
            hgexport(gcf,sprintf('%s',fstem{f_idx}), hgexport('factorystyle'),'Format','png');
    
        end;
end;

return;