close all; clear all;

subjects_dir='/Users/fhlin/workspace/subjects';
subject='fhlin';

file_surf={
{'outer_skin.surf'},
{'outer_skull.surf'},
{'inner_skull.surf'},
{'../surf/lh.orig'};
{'../surf/rh.orig'},
};

output_file_surf={
'fhlin_skin.mat'
'fhlin_skull.mat'
'fhlin_csf.mat'
'fhlin_gm_lh.mat'
'fhlin_gm_rh.mat'
};

for f_idx=1:length(file_surf)
    vv=[];ff=[];
    for i_idx=1:length(file_surf{f_idx})
        if(length(file_surf{f_idx})==1)
            [vv,ff]=read_surf(sprintf('%s/%s/bem/%s',subjects_dir,subject, file_surf{f_idx}{i_idx}));
        else
            [vv_tmp,ff_tmp]=read_surf(sprintf('%s/%s/%s',subjects_dir,subject, file_surf{f_idx}{i_idx}));
            vv=cat(1,vv,vv_tmp);
            ff=cat(1,ff,ff_tmp);
        end;
    end;
    surf_vertex{f_idx}=vv;
    surf_face{f_idx}=ff;

    surf_face{f_idx}=surf_face{f_idx}+1;
    
    TR = triangulation(surf_face{f_idx},surf_vertex{f_idx});
    surf_center{f_idx}= incenter(TR);
    surf_norm{f_idx} = faceNormal(TR);

    %save files for e-field modeling
    P=surf_vertex{f_idx};
    t=surf_face{f_idx};
    normals=surf_norm{f_idx};
    save(sprintf('%s',output_file_surf{f_idx}),'P','t','normals');

    hold on;
    colors=get(gca,'colororder');
    h=patch('vertices',surf_vertex{f_idx},'faces',surf_face{f_idx},'edgecolor','none','facecolor',colors(f_idx,:),'facealpha',0.2); 
    %quiver3(surf_center{f_idx}(:,1),surf_center{f_idx}(:,2),surf_center{f_idx}(:,3),surf_norm{f_idx}(:,1),surf_norm{f_idx}(:,2),surf_norm{f_idx}(:,3),0.5,'color','r');

end;
view(-160,20);
lighting phong
camlight
axis off vis3d equal