close all; clear all;


file_meg_mat='./meg_fwd_prep_032819.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('erf_041019.mat');

data=erf_all(1).erf.*1e15; %fT
timeVec=erf_all(1).timeVec*1e3; %ms
flag_baseline_correct=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(flag_baseline_correct)
    base_idx=find(timeVec<0);
    data=bsxfun(@minus,data,mean(data(:,base_idx),2));
end;

topo_value=data(:,300);
topo_threshold=[50  200];

load(file_meg_mat);
elec.elecpos=points;
elec.label=points_label;
n_meeg=size(points,1)/2;






for v_idx=1:n_meeg
    
    d=sqrt(sum(abs(bsxfun(@minus, verts_osc, points(v_idx,:))).^2,2));
    [dummy,topo_vertex(v_idx)]=min(d);
    
%    nn=points(v_idx+n_meeg,:)-points(v_idx,:);
%    nn=nn./norm(nn)./50;
    
%    h=quiver3(points(v_idx,1),points(v_idx,2),points(v_idx,3),nn(1),nn(2),nn(3));
%    set(h,'maxheadsize',1,'autoscale','off')

end;
figure(1);
etc_render_topo('vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_stc',data,'topo_vertex',topo_vertex,'topo_threshold',topo_threshold);
%etc_render_topo('vol_vertex',verts_osc,'vol_face',faces_osc-1);
hold on;

for v_idx=1:n_meeg
    
%    figure(1);
    
    h1=plot3(verts_osc(topo_vertex(v_idx),1),verts_osc(topo_vertex(v_idx),2),verts_osc(topo_vertex(v_idx),3),'r.');

    h2=plot3(points(v_idx,1),points(v_idx,2),points(v_idx,3),'g.');
    
%    figure(2);
%    v_idx
%    plot(data(v_idx,:))
    
    
%    if(v_idx>96)
%        keyboard;
%    end;
    
    delete(h1);
    delete(h2);

end;

topo_vertex=topo_vertex-1;

axis equal off vis3d;

hold on;

for v_idx=1:n_meeg
    %h=plot3(points(v_idx,1),points(v_idx,2),points(v_idx,3),'r.');
    h=plot3(verts_osc(topo_vertex(v_idx),1),verts_osc(topo_vertex(v_idx),2),verts_osc(topo_vertex(v_idx),3),'r.');
    
end;
return;





