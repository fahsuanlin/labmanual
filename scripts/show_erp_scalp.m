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
ERP=tmp./size(erp,1);
timeVec=erp{1,3}.timeVec;


%load('/Users/fhlin_admin/workspace/eegmri/180322_PYW/analysis/bem.mat');
%verts_osc_electrode_idx([32 33 34])=[];
load('eeg_fwd_prep_042319.mat');
for v_idx=1:length(points)
    dist=verts_osc-repmat(points(v_idx,:),[size(verts_osc,1) 1]);
    dist=sqrt(sum(dist.^2,2));
    [dummy,verts_osc_electrode_idx(v_idx)]=min(dist);
end;
etc_render_topo('vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_vertex',verts_osc_electrode_idx-1,'topo_stc',ERP,'topo_smooth',20,'topo_threshold',[1 2],'topo_stc_timevec',timeVec,'topo_stc_timevec_unit','s');
view(0,0);


return;

