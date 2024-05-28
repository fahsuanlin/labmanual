clear all; close all

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/eegmri/subjects/'); %for MAC/Linux

file_surface={
    '/Users/fhlin/workspace/eegmri/subjects/180411_PYY/bem/outer_skin.surf';
    '/Users/fhlin/workspace/eegmri/subjects/180411_PYY/bem/outer_skull.surf';
    '/Users/fhlin/workspace/eegmri/subjects/180411_PYY/bem/inner_skull.surf';
};

R=0.1; %10% of triangulation faces
colors=colororder;

for f_idx=1:length(file_surface)
    fprintf('reading [%s]...\n',file_surface{f_idx});
    [vertex_coords{f_idx}, faces{f_idx}] = read_surf(file_surface{f_idx});
    faces{f_idx}=faces{f_idx}+1; %1-based faces

    fprintf('decimating [%s] for %2.2f triangulation...\n',file_surface{f_idx}, R);
    [faces_reduced{f_idx}, vertex_coords_reduced{f_idx}] = reducepatch(faces{f_idx}, vertex_coords{f_idx}, R); 


    [pp,ff,ee]=fileparts(file_surface{f_idx});
    fn=sprintf('%s/%s_d10%s',pp,ff,ee);
    fprintf('saving decimated surface [%s] ...\n',fn);
    write_surf(fn, vertex_coords_reduced{f_idx}, faces_reduced{f_idx}); %0-based faces conversion automatically!!

    color_idx=mod(f_idx-1,size(colors,1))+1;
    figure(1); hold on;
    h=patch('faces',faces{f_idx},'vertices',vertex_coords{f_idx},'edgecolor','none','facealpha',0.2); axis off image vis3d; view(-130, 45);
    set(h,'facecolor',colors(color_idx,:));

    figure(2); hold on;
    h=patch('faces',faces_reduced{f_idx},'vertices',vertex_coords_reduced{f_idx},'edgecolor','none','facealpha',0.2); axis off image vis3d;  view(-130, 45);
    set(h,'facecolor',colors(color_idx,:));
    
end;