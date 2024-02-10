close all; clear all;

load field_e_fhlin.mat	

t0  = t(Indicator==objectnumber, :);

clear xx yy zz
for idx=1:3
    xx(:,idx)=P(t0(:,idx),1);
    yy(:,idx)=P(t0(:,idx),2);
    zz(:,idx)=P(t0(:,idx),3);
end;

coords_tms(:,1)=mean(xx,2);
coords_tms(:,2)=mean(yy,2);
coords_tms(:,3)=mean(zz,2);
    
ff=figure; set(ff,'visible','off');
setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri_wm/subjects');
subject='s001';
etc_render_fsbrain('subject',subject);
global etc_render_fsbrain
close(ff);
coords=etc_render_fsbrain.vertex_coords./1e3;

knn_idx=knnsearch(coords_tms,coords,'K',1);

%temp1=etc_threshold(Etotal,0.9999);
temp1=Etotal;
val=temp1(knn_idx);

clear etc_render_fsbrain;
figure;
etc_render_fsbrain('subject',subject,'overlay_value',val,'overlay_vertex',[0:length(val)-1],'overlay_threshold',[50 500],'overlay_exclude_fstem','exclude-lh.label')
view(90,0)


%%% show TMS coil
load('output_coil_data_fhlin.mat');
p1 = patch('vertices', Coil.P.*1e3, 'faces', Coil.t);
p1.FaceColor = [0.72 0.45 0.2];
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.20;

pointsline=pointsline.*1e3;

plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

