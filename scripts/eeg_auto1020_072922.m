clear all; close all

surf_outer_skin='/Users/fhlin/workspace/eegmri_emotion/subjects/s015/surf/lh.seghead.tri';

file_sensor_landmarks='sensor_landmarks.mat';
% a variable named 'sensor' with a field named 'coords' 
%   sensor.coords(1,:): coordinates for LPA
%   sensor.coords(2,:): coordinates for Nasion
%   sensor.coords(3,:): coordinates for RPA
%   sensor.coords(4,:): coordinates for Inion

output_fstem='eeg_auto1020_077922';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[verts_osc, faces_osc] = inverse_read_tri(surf_outer_skin);
verts_osc=verts_osc.*1e-3;


load(file_sensor_landmarks);

Inion=sensor.coords(4,:);
Nasion=sensor.coords(2,:);
LPA=sensor.coords(1,:);
RPA=sensor.coords(3,:);



figure;
etc_render_topo('alpha',0.8,'vol_vertex',verts_osc,'vol_face',faces_osc-1); 
hold on;

verts_osc_avg=mean(verts_osc,1);

maskNI = mask_NasionInion(verts_osc-verts_osc_avg, Nasion-verts_osc_avg, Inion-verts_osc_avg, .02,verts_osc_avg);  % 20 mm slab

[~,iCz] = max(verts_osc(:,3))
Cz_tmp = verts_osc(iCz,:);
arc_NI = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, Nasion-verts_osc_avg, Cz_tmp-verts_osc_avg, Inion-verts_osc_avg, maskNI);
arc_NI = arc_NI+verts_osc_avg;


%plot3(verts_osc(find(maskNI),1), verts_osc(find(maskNI),2), verts_osc(find(maskNI),3),'b.');
plot3(arc_NI(:,1), arc_NI(:,2), arc_NI(:,3),'r.');
h=plot3(LPA(1),LPA(2),LPA(3),'g.'); set(h,'markersize',24);
h=plot3(RPA(1),RPA(2),RPA(3),'g.'); set(h,'markersize',24);
h=plot3(Inion(1),Inion(2),Inion(3),'g.'); set(h,'markersize',24);
h=plot3(Nasion(1),Nasion(2),Nasion(3),'g.'); set(h,'markersize',24);

elect.FPz = arc_interpolate(arc_NI, 10);
elect.AFz = arc_interpolate(arc_NI, 20);
elect.Fz  = arc_interpolate(arc_NI, 30);
elect.FCz  = arc_interpolate(arc_NI, 40);
elect.CPz  = arc_interpolate(arc_NI, 60);
elect.Pz  = arc_interpolate(arc_NI, 70);
elect.POz  = arc_interpolate(arc_NI, 80);
elect.Oz  = arc_interpolate(arc_NI, 90);
Cz  = arc_interpolate(arc_NI, 50);


maskLR = mask_LPARPA(verts_osc-verts_osc_avg, LPA-verts_osc_avg, RPA-verts_osc_avg, .02); % 20 mm slab
arc_LR = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, LPA-verts_osc_avg, Cz-verts_osc_avg, RPA-verts_osc_avg, maskLR);
arc_LR = arc_LR+verts_osc_avg;

%plot3(verts_osc(find(maskLR),1), verts_osc(find(maskLR),2), verts_osc(find(maskLR),3),'b.');
plot3(arc_LR(:,1), arc_LR(:,2), arc_LR(:,3),'r.');


elect.T7 = arc_interpolate(arc_LR, 10);  
elect.C5 = arc_interpolate(arc_LR, 20);  % exactly Cz plane
elect.C3 = arc_interpolate(arc_LR, 30);
elect.C1= arc_interpolate(arc_LR, 40);
elect.Cz = arc_interpolate(arc_LR, 50);
elect.T8 = arc_interpolate(arc_LR, 90);  
elect.C6 = arc_interpolate(arc_LR, 80);  % exactly Cz plane
elect.C4 = arc_interpolate(arc_LR, 70);
elect.C2= arc_interpolate(arc_LR, 60);
elect.Cz = arc_interpolate(arc_LR, 50);


arc_RR = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.Oz-verts_osc_avg, elect.T8-verts_osc_avg, elect.FPz-verts_osc_avg);
arc_RR=arc_RR+verts_osc_avg;
plot3(arc_RR(:,1), arc_RR(:,2), arc_RR(:,3),'r.');

arc_LL = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.Oz-verts_osc_avg, elect.T7-verts_osc_avg, elect.FPz-verts_osc_avg);
arc_LL=arc_LL+verts_osc_avg;
plot3(arc_LL(:,1), arc_LL(:,2), arc_LL(:,3),'r.');

elect.FT7= arc_interpolate(arc_LL, 60);
elect.FT8= arc_interpolate(arc_RR, 60);
elect.F7= arc_interpolate(arc_LL, 70);
elect.F8= arc_interpolate(arc_RR, 70);
elect.AF7= arc_interpolate(arc_LL, 80);
elect.AF8= arc_interpolate(arc_RR, 80);
elect.Fp1= arc_interpolate(arc_LL, 90);
elect.Fp2= arc_interpolate(arc_RR, 90);

elect.TP7= arc_interpolate(arc_LL, 40);
elect.TP8= arc_interpolate(arc_RR, 40);
elect.P7= arc_interpolate(arc_LL, 30);
elect.P8= arc_interpolate(arc_RR, 30);
elect.PO7= arc_interpolate(arc_LL, 20);
elect.PO8= arc_interpolate(arc_RR, 20);
elect.O1= arc_interpolate(arc_LL, 10);
elect.O2= arc_interpolate(arc_RR, 10);

arc_Fmid = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.F7-verts_osc_avg, elect.Fz-verts_osc_avg, elect.F8-verts_osc_avg);
arc_Fmid=arc_Fmid+verts_osc_avg;
plot3(arc_Fmid(:,1), arc_Fmid(:,2), arc_Fmid(:,3),'r.');

arc_Fant = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.AF7-verts_osc_avg, elect.AFz-verts_osc_avg, elect.AF8-verts_osc_avg);
arc_Fant=arc_Fant+verts_osc_avg;
plot3(arc_Fant(:,1), arc_Fant(:,2), arc_Fant(:,3),'r.');

arc_Fpos = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.FT7-verts_osc_avg, elect.FCz-verts_osc_avg, elect.FT8-verts_osc_avg);
arc_Fpos=arc_Fpos+verts_osc_avg;
plot3(arc_Fpos(:,1), arc_Fpos(:,2), arc_Fpos(:,3),'r.');

arc_Pmid = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.P7-verts_osc_avg, elect.Pz-verts_osc_avg, elect.P8-verts_osc_avg);
arc_Pmid=arc_Pmid+verts_osc_avg;
plot3(arc_Pmid(:,1), arc_Pmid(:,2), arc_Pmid(:,3),'r.');

arc_Pant = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.TP7-verts_osc_avg, elect.CPz-verts_osc_avg, elect.TP8-verts_osc_avg);
arc_Pant=arc_Pant+verts_osc_avg;
plot3(arc_Pant(:,1), arc_Pant(:,2), arc_Pant(:,3),'r.');

arc_Ppos = geodesic_arc_through_Cz(verts_osc-verts_osc_avg, faces_osc, elect.PO7-verts_osc_avg, elect.POz-verts_osc_avg, elect.PO8-verts_osc_avg);
arc_Ppos=arc_Ppos+verts_osc_avg;
plot3(arc_Ppos(:,1), arc_Ppos(:,2), arc_Ppos(:,3),'r.');

elect.AF3= arc_interpolate(arc_Fant, 25);
elect.AF4= arc_interpolate(arc_Fant, 75);

elect.F5= arc_interpolate(arc_Fmid, 12.5);
elect.F3= arc_interpolate(arc_Fmid, 25);
elect.F1= arc_interpolate(arc_Fmid, 37.5);
elect.F2= arc_interpolate(arc_Fmid, 62.5);
elect.F4= arc_interpolate(arc_Fmid, 75);
elect.F6= arc_interpolate(arc_Fmid, 87.5);

elect.FC5= arc_interpolate(arc_Fpos, 12.5);
elect.FC3= arc_interpolate(arc_Fpos, 25);
elect.FC1= arc_interpolate(arc_Fpos, 37.5);
elect.FC2= arc_interpolate(arc_Fpos, 62.5);
elect.FC4= arc_interpolate(arc_Fpos, 75);
elect.FC6= arc_interpolate(arc_Fpos, 87.5);

elect.PO3= arc_interpolate(arc_Ppos, 25);
elect.PO4= arc_interpolate(arc_Ppos, 75);

elect.P5= arc_interpolate(arc_Pmid, 12.5);
elect.P3= arc_interpolate(arc_Pmid, 25);
elect.P1= arc_interpolate(arc_Pmid, 37.5);
elect.P2= arc_interpolate(arc_Pmid, 62.5);
elect.P4= arc_interpolate(arc_Pmid, 75);
elect.P6= arc_interpolate(arc_Pmid, 87.5);

elect.CP5= arc_interpolate(arc_Pant, 12.5);
elect.CP3= arc_interpolate(arc_Pant, 25);
elect.CP1= arc_interpolate(arc_Pant, 37.5);
elect.CP2= arc_interpolate(arc_Pant, 62.5);
elect.CP4= arc_interpolate(arc_Pant, 75);
elect.CP6= arc_interpolate(arc_Pant, 87.5);

labels = fieldnames(elect);
for k = 1:length(labels)
    p(k,:) = elect.(labels{k});
end
labels{end+1}='LPA';
labels{end+1}='RPA';
labels{end+1}='Nasion';
labels{end+1}='Inion';
p(end+1,:)=LPA;
p(end+1,:)=RPA;
p(end+1,:)=Nasion;
p(end+1,:)=Inion;

figure;
etc_render_topo('alpha',0.8,'vol_vertex',verts_osc,'vol_face',faces_osc-1,'topo_aux_point_coords',p,'topo_aux_point_name',labels);
view(-60,20);

save(sprintf('%s.mat',output_fstem),'p','labels');

return;













function p = arc_interpolate(arc, percent)
    % arc: Nx3
    d = vecnorm(diff(arc),2,2);
    cum = [0; cumsum(d)];
    total = cum(end);

    dist = percent/100 * total;
    idx = find(cum >= dist,1);

    % linear interpolate
    if idx==1
        p = arc(1,:);
    else
        t = (dist - cum(idx-1)) / (cum(idx) - cum(idx-1));
        p = arc(idx-1,:) + t*(arc(idx,:) - arc(idx-1,:));
    end
end



function mask = mask_NasionInion(verts_centered, Nasion, Inion, slab_halfwidth,varargin)

    if(length(varargin)>0)
        verts_avg=varargin{1};
    end;

    % verts_centered: Nx3 mesh centered at origin
    % slab_halfwidth: e.g., 15–25 mm
    
    X = verts_centered(:,1);
    
    % keep vertices where x≈0 (slab) and z>0 (upper hemisphere)
    %mask = abs(X) < slab_halfwidth & Z > 0;
    idx_band = find(abs(X) < slab_halfwidth);

%    plot3(verts_centered((idx_band),1)+verts_avg(1), verts_centered((idx_band),2)+verts_avg(2), verts_centered((idx_band),3)+verts_avg(3),'c.');

    Vband = verts_centered(idx_band,:);

    dN = vecnorm(Vband - Nasion, 2, 2);
    dI = vecnorm(Vband - Inion,  2, 2);

    isNgroup = dN < dI;     % logical vector
    isIgroup = ~isNgroup;

    Znas = Nasion(3);
    Zin  = Inion(3);

    keep_N = isNgroup & (Vband(:,3) >= Znas);
    keep_I = isIgroup & (Vband(:,3) >= Zin);

    keep_band = keep_N | keep_I;

    % Convert 'keep_band' back to full-vertex logical mask
    mask = false(size(verts_centered,1),1);
    mask(idx_band(keep_band)) = true;
    

%    plot3(verts_centered(idx_band(find(keep_N)),1)+verts_avg(1), verts_centered(idx_band(find(keep_N)),2)+verts_avg(2), verts_centered(idx_band(find(keep_N)),3)+verts_avg(3),'r.');
%    plot3(verts_centered((mask),1)+verts_avg(1), verts_centered((mask),2)+verts_avg(2), verts_centered((mask),3)+verts_avg(3),'b.');

end

function mask = mask_LPARPA(verts_centered, LPA, RPA, slab_halfwidth)
    Y = verts_centered(:,2);
    
    % x-z plane slab around y=0, only upper scalp
    %mask = abs(Y) < slab_halfwidth & Z > 0;
    idx_band =find(abs(Y) < slab_halfwidth);
    Vband = verts_centered(idx_band,:);

    dL = vecnorm(Vband - LPA, 2, 2);
    dR = vecnorm(Vband - RPA,  2, 2);

    isLgroup = dL < dR;     % logical vector
    isRgroup = ~isLgroup;

    ZLPA = LPA(3);
    ZRPA  = RPA(3);

    keep_L = isLgroup & (Vband(:,3) >= ZLPA);
    keep_R = isRgroup & (Vband(:,3) >= ZRPA);

    keep_band = keep_L | keep_R;

    % Convert 'keep_band' back to full-vertex logical mask
    mask = false(size(verts_centered,1),1);
    mask(idx_band(keep_band)) = true;
end

function arc = geodesic_arc_through_Cz(verts_c, faces, p_start, p_Cz, p_end, varargin)
    if(length(varargin)>0)
        mask=varargin{1};
    else
        mask=ones(size(verts_c,1),1);
    end;

    % segment 1: start -> Cz
    arc1 = eeg_constrained_geodesic(verts_c, faces, p_start, p_Cz, mask);

    % segment 2: Cz -> end
    arc2 = eeg_constrained_geodesic(verts_c, faces, p_Cz, p_end, mask);

    % concatenate, drop duplicate Cz vertex at the junction
    arc = [arc1; arc2(2:end,:)];
end