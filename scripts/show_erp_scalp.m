close all; clear all;

load erp_avg_outside.mat;
%read the 3rd condition; this is subject to change
erp_condition_index=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=[];
for i=1:size(erp_avg,1)
     if(isempty(tmp))
        tmp=erp_avg{i,erp_condition_index}.erp;
    else
        tmp=tmp+erp_avg{i,erp_condition_index}.erp;
    end;
end;
vep=tmp./size(erp_avg,1);
timeVec=erp_avg{1,erp_condition_index}.timeVec;

%31-channel EEG topography definition and a head model for Brain Products
load topology_31ch_default.mat;
etc_render_topo('vol_vertex',vertex,'vol_face',face-1,'topo_vertex',electrode_idx-1,'topo_stc',vep,'topo_smooth',10,'topo_threshold',[5 10],'topo_stc_timevec',timeVec,'topo_stc_timevec_unit','s');


return;

