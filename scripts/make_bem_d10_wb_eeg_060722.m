clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri/subjects/'); %for MAC/Linux

%source space
file_source_fif='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/180322_PYW-5-src.fif';
file_source_wholebrain='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/mri/aseg.mgz';
file_source_wholebrain_orig='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/mri/orig.mgz';

wholebrain_index={
    [6 8 9 10 11 12 13 16 17 18 19 20 26 27 ];
    [45 47 48 49 50 51 52 53 54 55 56 58 59];
    };
%the following indices were from FreeSurfer look-up table 
% 6   Left-Cerebellum-Exterior                0   148 0   0
% 8   Left-Cerebellum-Cortex                  230 148 34  0
% 9   Left-Thalamus                           0   118 14  0
% 10  Left-Thalamus-Proper*                   0   118 14  0
% 11  Left-Caudate                            122 186 220 0
% 12  Left-Putamen                            236 13  176 0
% 13  Left-Pallidum                           12  48  255 0
% 16  Brain-Stem                              119 159 176 0
% 17  Left-Hippocampus                        220 216 20  0
% 18  Left-Amygdala                           103 255 255 0
% 19  Left-Insula                             80  196 98  0
% 20  Left-Operculum                          60  58  210 0
% 26  Left-Accumbens-area                     255 165 0   0
% 27  Left-Substancia-Nigra                   0   255 127 0
% 45  Right-Cerebellum-Exterior               0   148 0   0
% 47  Right-Cerebellum-Cortex                 230 148 34  0
% 48  Right-Thalamus                          0   118 14  0
% 49  Right-Thalamus-Proper*                  0   118 14  0
% 50  Right-Caudate                           122 186 220 0
% 51  Right-Putamen                           236 13  176 0
% 52  Right-Pallidum                          13  48  255 0
% 53  Right-Hippocampus                       220 216 20  0
% 54  Right-Amygdala                          103 255 255 0
% 55  Right-Insula                            80  196 98  0
% 56  Right-Operculum                         60  58  210 0
% 58  Right-Accumbens-area                    255 165 0   0
% 59  Right-Substancia-Nigra                  0   255 127 0


surf_outer_skin='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/outer_skin_d10.surf';
surf_outer_skull='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/outer_skull_d10.surf';
surf_inner_skull='/Users/fhlin/workspace/eegmri/subjects/180322_PYW/bem/inner_skull_d10.surf';

%sesnor info
file_eeg_mat='./eeg_fwd_prep_060722.mat';

select_channel={ 'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'};


output_stem='180322_PYW_d10_wb_eeg_060722';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%LOADING THE SOURCE SPACE%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[src] = mne_read_source_spaces(file_source_fif);

%whole brain source space
if(~isempty(file_source_wholebrain_orig))
    src_wb_orig=MRIread(file_source_wholebrain_orig);
end;
if(~isempty(file_source_wholebrain))
    src_wb=MRIread(file_source_wholebrain);
    src_wb_vol=src_wb.vol;
    
    for hemi_idx=1:2
        src_wb_idx{hemi_idx}=[];
        wb_index{hemi_idx}=[];
        for idx=1:length(wholebrain_index{hemi_idx})
            tmp=find(src_wb_vol(:)==wholebrain_index{hemi_idx}(idx));
            src_wb_idx{hemi_idx}=cat(1,src_wb_idx{hemi_idx},tmp(:));
            wb_index{hemi_idx}=cat(1,wb_index{hemi_idx},ones(size(tmp(:))).*wholebrain_index{hemi_idx}(idx));
        end;
        
        for idx=1:length(src_wb_idx{hemi_idx})
            [cc,rr,ss]=ind2sub(size(src_wb_orig.vol),src_wb_idx{hemi_idx}(idx));
            tmp=src_wb.tkrvox2ras*[rr cc ss 1]';
            wb_coord{hemi_idx}(idx,:)=tmp(1:3)'; %this is the coordinate of whole-brain source in the surface coordinates
        end;
        
        [cc,rr,ss]=ind2sub(size(src_wb_orig.vol),src_wb_idx{hemi_idx}(:));
        cc_min=min(cc(:));
        cc_max=max(cc(:));
        rr_min=min(rr(:));
        rr_max=max(rr(:));
        ss_min=min(ss(:));
        ss_max=max(ss(:));
        cc_idx=ismember(cc(:),[cc_min:2:cc_max]);
        rr_idx=ismember(rr(:),[rr_min:2:rr_max]);
        ss_idx=ismember(ss(:),[ss_min:2:ss_max]);
        crs_idx=(cc_idx&rr_idx&ss_idx);
        src_wb_idx{hemi_idx}=src_wb_idx{hemi_idx}(crs_idx);
        wb_coord{hemi_idx}=wb_coord{hemi_idx}(crs_idx,:);
%         rr=rr(find(crs_idx));
%         cc=cc(find(crs_idx));
%         ss=ss(find(crs_idx));
    end;
end;


[verts_osc, faces_osc] = mne_read_surface(surf_outer_skin);
[verts_osk, faces_osk] = mne_read_surface(surf_outer_skull);
[verts_isk, faces_isk] = mne_read_surface(surf_inner_skull);
verts_osc=verts_osc.*1e3; %in to mm
verts_osk=verts_osk.*1e3; %in to mm
verts_isk=verts_isk.*1e3; %in to mm

%verts_isk(:,3)=verts_isk(:,3)+5;

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
h_isk=patch('Vertices',verts_osc,'Faces',faces_osc,'Facealpha',0.1,'FaceColor',[1 0 1],'EdgeColor','none');
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


tri = delaunayn([verts_isk(:,1) verts_isk(:,2) verts_isk(:,3)]); % Generate delaunay triangulization
        

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
    v=[verts_osc(v_idx,1),verts_osc(v_idx,2),verts_osc(v_idx,3)];
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
%scalp
verts_osk_norm=zeros(size(verts_osk));
for v_idx=1:size(verts_osk,1)
    tmp=[verts_osk(v_idx,1),verts_osc(v_idx,2),verts_osc(v_idx,3)];
    verts_osk_norm(v_idx,:)=tmp(:).'/norm(tmp);
end;




Index=find(contains(points_label,select_channel));
if(length(Index)<=length(select_channel)) %all electrodes were found on topology
    for ii=1:length(select_channel)
        for idx=1:length(points_label)
            if(strcmp(points_label{idx},select_channel{ii}))
                Index(ii)=idx;
                electrode_data_idx(idx)=ii;
            end;
        end;
    end;
    if(isempty(Index))
        fprintf('cannot find corresponding channels!\nerror in loading the topology!\n');
        
        return;
    end;
end;
                                            
                                            
%IDX = knnsearch(X,Y) 
for v_idx=1:n_eeg
    
    fprintf(fp,'%f %f %f\n',points(Index(v_idx),1).*1e3, points(Index(v_idx),2).*1e3,points(Index(v_idx),3).*1e3);
  
end;
hold on;
hh=plot3(points(Index(1:n_eeg),1).*1e3, points(Index(1:n_eeg),2).*1e3,points(Index(1:n_eeg),3).*1e3,'r.');
set(hh,'markersize',12);
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
   
   for v_idx=1:size(wb_coord{s_idx},1)
       %if(wb_coord{s_idx}(v_idx,3)>-66)
           fprintf(fp,'%f %f %f %f %f %f\n',wb_coord{s_idx}(v_idx,1),wb_coord{s_idx}(v_idx,2),wb_coord{s_idx}(v_idx,3),1,0,0);
           fprintf(fp,'%f %f %f %f %f %f\n',wb_coord{s_idx}(v_idx,1),wb_coord{s_idx}(v_idx,2),wb_coord{s_idx}(v_idx,3),0,1,0);
           fprintf(fp,'%f %f %f %f %f %f\n',wb_coord{s_idx}(v_idx,1),wb_coord{s_idx}(v_idx,2),wb_coord{s_idx}(v_idx,3),0,0,1);
       %end;
   end;
   
   hold on;
   %plot3(src(s_idx).rr(src(s_idx).vertno(:),1).*1e3, src(s_idx).rr(src(s_idx).vertno(:),2).*1e3, src(s_idx).rr(src(s_idx).vertno(:),3).*1e3,'r.');
   cort{s_idx}(:,1)=src(s_idx).rr(src(s_idx).vertno(:),1).*1e3;
   cort{s_idx}(:,2)=src(s_idx).rr(src(s_idx).vertno(:),2).*1e3;
   cort{s_idx}(:,3)=src(s_idx).rr(src(s_idx).vertno(:),3).*1e3;
   
   %idx=find(wb_coord{s_idx}(:,3)>-66);
   %plot3(wb_coord{s_idx}(idx,1),wb_coord{s_idx}(idx,2),wb_coord{s_idx}(idx,3),'g.');
   plot3(wb_coord{s_idx}(:,1),wb_coord{s_idx}(:,2),wb_coord{s_idx}(:,3),'g.');
   fclose(fp);
end;

%save bem_d10_wb_091019.mat verts_osc verts_osk verts_isk faces_osc faces_isk faces_osk wb_coord cort 

return;
