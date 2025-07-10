close all; clear all;

file_entorhinal={
    'entorhinal_dfconn_native_vol_aseg_060525_gavg_entorhinal-anat.nii';
    };

file_hippo={
    'hippo_dfconn_native_vol_aseg_060525_gavg_hippo_left-anat.nii';
    'hippo_dfconn_native_vol_aseg_060525_gavg_hippo_right-anat.nii';
    };


for f_idx=1:length(file_entorhinal)
    tmp=MRIread(file_entorhinal{f_idx});
    if(f_idx==1)
        ent=tmp.vol;
    else
        ent=ent+tmp.vol;
    end;
end;
ent=reshape(ent,[size(ent,1)*size(ent,2)*size(ent,3),size(ent,4)]);


for f_idx=1:length(file_hippo)
    tmp=MRIread(file_hippo{f_idx});
    if(f_idx==1)
        hippo=tmp.vol;
    else
        hippo=hippo+tmp.vol;
    end;
end;
hippo=reshape(hippo,[size(hippo,1)*size(hippo,2)*size(hippo,3),size(hippo,4)]);

[pval,tstat]=etc_ttest2(ent',hippo');
% 
% for v_idx=1:size(ent,1)
%      [H,P,CI,STATS] = ttest2(ent(v_idx,:),hippo(v_idx,:));
%      tstat(v_idx)=STATS.tstat;
% end;

tstat=reshape(tstat,[size(tmp.vol,1),size(tmp.vol,2),size(tmp.vol,3)]);
tmp.nframes=1;
tmp.vol=tstat;
MRIwrite(tmp,'compare_dfconn_hippo_entorhinal_060525-anat.nii');