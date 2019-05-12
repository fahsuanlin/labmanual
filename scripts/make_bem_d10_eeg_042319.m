clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/sinica_meg/subjects/'); %for MAC/Linux

%source space
file_source_fif='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/041019-5-src.fif';

surf_outer_skin='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/outer_skin_d10.surf';
surf_outer_skull='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/outer_skull_d10.surf';
surf_inner_skull='/Users/fhlin_admin/workspace/sinica_meg/subjects/041019/bem/inner_skull_d10.surf';

%sesnor info
file_eeg_mat='./eeg_fwd_prep_042319.mat';

output_stem='041019_d10';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[src] = mne_read_source_spaces(file_source_fif);


load(file_eeg_mat);
elec.elecpos=points;
elec.label=points_label;


[verts_osc, faces_osc] = mne_read_surface(surf_outer_skin);
[verts_osk, faces_osk] = mne_read_surface(surf_outer_skull);
[verts_isk, faces_isk] = mne_read_surface(surf_inner_skull);
verts_osc=verts_osc.*1e3; %in to mm
verts_osk=verts_osk.*1e3; %in to mm
verts_isk=verts_isk.*1e3; %in to mm



src_stem={'lh','rh'};
xx=[];
yy=[];
zz=[];
for s_idx=1:length(src)
    fp=fopen(sprintf('%s-%s.dip',output_stem,src_stem{s_idx}),'w');
    xx=cat(1,xx,src(s_idx).rr(:,1));
    yy=cat(1,yy,src(s_idx).rr(:,2));
    zz=cat(1,zz,src(s_idx).rr(:,3));
    fclose(fp);
end;
h=plot3(xx,yy,zz,'g.');

hold on;

points=[];
%etc_render_topo('vol_vertex',verts_osc,'vol_face',faces_osc-1,'alpha',0.3);
hold on;
h_osc=patch('Vertices',verts_osc,'Faces',faces_osc,'Facealpha',0.1,'FaceColor',[1 0 1],'EdgeColor','none');
h_isk=patch('Vertices',verts_isk,'Faces',faces_isk,'Facealpha',0.1,'FaceColor',[1 1 0],'EdgeColor','none');
h_osk=patch('Vertices',verts_osk,'Faces',faces_osk,'Facealpha',0.1,'FaceColor',[0 1 1],'EdgeColor','none');

axis vis3d equal 
lighting phong
camlight


%make .geom file
fprintf('making geometry file...\n');
fp=fopen(sprintf('%s_bem.geom',output_stem),'w');
fprintf(fp,'# Domain Description 1.1\n');
fprintf(fp,'Interfaces 3\n');
fprintf(fp,'Interface cortex: "cortex.tri"\n');
fprintf(fp,'Interface skull: "skull.tri"\n');
fprintf(fp,'Interface scalp: "scalp.tri"\n');
fprintf(fp,'Domains 4\n');
fprintf(fp,'Domain brain: -cortex\n');
fprintf(fp,'Domain skull: -skull +cortex\n');
fprintf(fp,'Domain scalp: -scalp +skull\n');
fprintf(fp,'Domain air: +scalp\n');
fclose(fp);

%make conductivity file
fprintf('making conductivity file...\n');
fp=fopen(sprintf('%s_bem.cond',output_stem),'w');
fprintf(fp,'# Properties Description 1.0 (Conductivities)\n');
fprintf(fp,'air 0.0\n');
fprintf(fp,'scalp 1\n');
fprintf(fp,'brain 1\n');
fprintf(fp,'skull 0.0125\n');
fclose(fp);

%make mesh files
fprintf('making mesh files...\n');

fprintf('making inner skull mesh file...\n');
fp=fopen(sprintf('cortex.tri'),'w');
fprintf(fp,'- %d\n',size(verts_isk,1));
for v_idx=1:size(verts_isk,1)
    v=[verts_isk(v_idx,1),verts_isk(v_idx,2),verts_isk(v_idx,3)];
    v=v./norm(v);
    fprintf(fp,'%f %f %f \n',verts_isk(v_idx,1),verts_isk(v_idx,2),verts_isk(v_idx,3),v(1),v(2),v(3));
end;
fprintf(fp,'- %d %d %d\n',size(faces_isk,1), size(faces_isk,1), size(faces_isk,1));
for f_idx=1:size(faces_isk,1)
    fprintf(fp,'%d %d %d \n',faces_isk(f_idx,1)-1,faces_isk(f_idx,2)-1,faces_isk(f_idx,3)-1);
end;
fclose(fp);

fprintf('making outer skull mesh file...\n');
fp=fopen(sprintf('skull.tri'),'w');
fprintf(fp,'- %d\n',size(verts_osk,1));
for v_idx=1:size(verts_osk,1)
    v=[verts_osk(v_idx,1),verts_osk(v_idx,2),verts_osk(v_idx,3)];
    v=v./norm(v);
    fprintf(fp,'%f %f %f \n',verts_osk(v_idx,1),verts_osk(v_idx,2),verts_osk(v_idx,3),v(1),v(2),v(3));
end;
fprintf(fp,'- %d %d %d\n',size(faces_osk,1), size(faces_osk,1), size(faces_osk,1));
for f_idx=1:size(faces_osk,1)
    fprintf(fp,'%d %d %d \n',faces_osk(f_idx,1)-1,faces_osk(f_idx,2)-1,faces_osk(f_idx,3)-1);
end;
fclose(fp);

fprintf('making outer skin mesh file...\n');
fp=fopen(sprintf('scalp.tri'),'w');
fprintf(fp,'- %d\n',size(verts_osc,1));
for v_idx=1:size(verts_osk,1)
    v=[verts_osk(v_idx,1),verts_osc(v_idx,2),verts_osc(v_idx,3)];
    v=v./norm(v);
    fprintf(fp,'%f %f %f \n',verts_osc(v_idx,1),verts_osc(v_idx,2),verts_osc(v_idx,3),v(1),v(2),v(3));
end;
fprintf(fp,'- %d %d %d\n',size(faces_osc,1), size(faces_osc,1), size(faces_osc,1));
for f_idx=1:size(faces_osc,1)
    fprintf(fp,'%d %d %d \n',faces_osc(f_idx,1)-1,faces_osc(f_idx,2)-1,faces_osc(f_idx,3)-1);
end;
fclose(fp);

%make EEG sensor file
fprintf('making EEG sensor file...\n');
load(file_eeg_mat); %including variable 'n_meg'

fp=fopen(sprintf('%s.eegsensors',output_stem),'w');
for v_idx=1:n_eeg
    
    fprintf(fp,'%f %f %f\n',points(v_idx,1).*1e3, points(v_idx,2).*1e3,points(v_idx,3).*1e3);
  
end;
hold on;
plot3(points(:,1).*1e3, points(:,2).*1e3,points(:,3).*1e3,'b.');
fclose(fp);

%make source files
fprintf('making source files...\n');
src_stem={'lh','rh'};
for s_idx=1:length(src)
    fp=fopen(sprintf('%s-%s.dip',output_stem,src_stem{s_idx}),'w');
    for v_idx=1:src(s_idx).nuse
        fprintf(fp,'%f %f %f %f %f %f\n',src(s_idx).rr(src(s_idx).vertno(v_idx),1).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),2).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),3).*1e3, 1,0,0); 
        fprintf(fp,'%f %f %f %f %f %f\n',src(s_idx).rr(src(s_idx).vertno(v_idx),1).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),2).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),3).*1e3, 0,1,0); 
        fprintf(fp,'%f %f %f %f %f %f\n',src(s_idx).rr(src(s_idx).vertno(v_idx),1).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),2).*1e3, src(s_idx).rr(src(s_idx).vertno(v_idx),3).*1e3, 0,0,1); 
    end;
    
    hold on;
    plot3(src(s_idx).rr(src(s_idx).vertno(:),1).*1e3, src(s_idx).rr(src(s_idx).vertno(:),2).*1e3, src(s_idx).rr(src(s_idx).vertno(:),3).*1e3,'r.');


    fclose(fp);
end;



return;
