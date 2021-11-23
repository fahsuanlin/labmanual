close all; clear all;

%subject
subject='s031';

%electrode coodinates for *pre-surgery* MRI
file_electrode='electrode_101721_223845_s031.mat';

%electrode data
file_erf='erf_061721.mat';

%labels
file_label={
    'aud_HCPMMP1-lh_s031.label';
    'aud_a1_HCPMMP1-lh_s031.label';
    'aud_belt_HCPMMP1-lh_s031.label';
    'aud_HCPMMP1-rh_s031.label';
    'aud_a1_HCPMMP1-rh_s031.label';
    'aud_belt_HCPMMP1-rh_s031.label';
    'lh.G_temp_sup-G_T_transv_s031.label';
    'lh.G_temp_sup-Lateral_div1_s031.label';
    'lh.G_temp_sup-Lateral_div2_s031.label';
    'lh.planum_temporale_s031.label';
    'lh.S_intrapariet_and_P_trans_s031.label';
    'lh.S_occipital_ant_s031.label';
    'lh.S_temporal_sup_div1_s031.label';
    'lh.S_temporal_sup_div2_s031.label';
    'lh.S_temporal_sup_div3_s031.label';
    'lh.S_temporal_transverse_s031.label';
    'rh.G_temp_sup-G_T_transv_s031.label';
    'rh.G_temp_sup-Lateral_div1_s031.label';
    'rh.G_temp_sup-Lateral_div2_s031.label';
    'rh.planum_temporale_s031.label';
    'rh.S_intrapariet_and_P_trans_s031.label';
    'rh.S_occipital_ant_s031.label';
    'rh.S_temporal_sup_div1_s031.label';
    'rh.S_temporal_sup_div2_s031.label';
    'rh.S_temporal_sup_div3_s031.label';
    'rh.S_temporal_transverse_s031.label';
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load surface coordinates (orig.mgz)
subjects_dir=getenv('SUBJECTS_DIR');
file_orig_surf=sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,'lh','orig');
[orig_vertex_coords_lh, orig_faces_lh] = read_surf(file_orig_surf);
file_orig_surf=sprintf('%s/%s/surf/%s.%s',subjects_dir,subject,'rh','orig');
[orig_vertex_coords_rh, orig_faces_rh] = read_surf(file_orig_surf);

%%% load labels
for label_idx=1:length(file_label)
        [dummy,fstem]=fileparts(file_label{label_idx});
        fprintf('loading <%s>...',file_label{label_idx});
        [roi(label_idx).idx]=inverse_read_label(file_label{label_idx},'flag_display',0);
        roi(label_idx).name=fstem;
        if(isempty(findstr(fstem,'lh.'))&&isempty(findstr(fstem,'-lh')))
                fprintf('[RH]...\n');
                roi(label_idx).hemi='rh';
        else
                fprintf('[LH]...\n');
                roi(label_idx).hemi='lh';
        end;
end;

%%% load electrode info and coordinates
load(file_electrode);


%%% load electrode data
if(~isempty(file_erf))
    load(file_erf);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for label_idx=1:length(file_label)
    
    switch(roi(label_idx).hemi)
        case 'lh'
            label_coords=orig_vertex_coords_lh(roi(label_idx).idx,:);
        case 'rh'
            label_coords=orig_vertex_coords_rh(roi(label_idx).idx,:);
    end;
    label_coords_com=mean(label_coords,1); %center of mass
    
    %find electrode contacts closest to the selected label
    if(~isempty(electrode))
        
        max_contact=0;
        for e_idx=1:length(electrode)
            if(electrode(e_idx).n_contact>max_contact)
                max_contact=electrode(e_idx).n_contact;
            end;
        end;
        electrode_dist_min=ones(length(electrode),max_contact).*nan;
        electrode_dist_avg=ones(length(electrode),max_contact).*nan;
        
        for e_idx=1:length(electrode)
            for c_idx=1:electrode(e_idx).n_contact
                
                surface_coord=electrode(e_idx).coord(c_idx,:);
                
                tmp=label_coords-repmat(surface_coord(:)',[size(label_coords,1),1]);
                tmp=sqrt(sum(tmp.^2,2));
                
                electrode_dist_min(e_idx,c_idx)=min(tmp);
                electrode_dist_avg(e_idx,c_idx)=mean(tmp);                
            end;
        end;
        
        electrode_dist_min_com=ones(length(electrode),max_contact).*nan;
        
        for e_idx=1:length(electrode)
            for c_idx=1:electrode(e_idx).n_contact
                
                surface_coord=electrode(e_idx).coord(c_idx,:);
                
                tmp=surface_coord-label_coords_com;
                tmp=sqrt(sum(tmp.^2,2));
                
                electrode_dist_min_com(e_idx,c_idx)=tmp;
            end;
        end;
        
        
        [dummy,min_idx]=sort(electrode_dist_min(:));
        fprintf('<<%s>>\n',roi(label_idx).name);
        for ii=1:3 %show the nearest three contacts
            [ee,cc]=ind2sub(size(electrode_dist_min),min_idx(ii));
            fprintf('\tmin dist::<%s_%02d> %2.2f (mm) (%1.1f %1.1f %1.1f)\n',electrode(ee).name,cc,dummy(ii),electrode(ee).coord(cc,1),electrode(ee).coord(cc,2),electrode(ee).coord(cc,3));
        end;
        flag_found=0;
        electrode_idx=1;
        while(~flag_found)
            [ee,cc]=ind2sub(size(electrode_dist_min),min_idx(electrode_idx));
            roi(label_idx).electrode_min_dist_electrode_name=electrode(ee).name;
            roi(label_idx).electrode_min_dist_electrode_contact=cc;
            roi(label_idx).electrode_min_dist=dummy(1);
            target_electrode_contact=sprintf('%s%d',electrode(ee).name,cc);
            IndexC = strfind(erf_all(1).name,target_electrode_contact);
            Index = find(not(cellfun('isempty',IndexC)));
            if(~isempty(Index))
                flag_found=1;
                for cond_idx=1:length(erf_all)
                    roi(label_idx).erf_electrode_min_dist(cond_idx).data=squeeze(erf_all(cond_idx).erf_raw(Index,:,:));
                end;
            else
                electrode_idx=electrode_idx+1;
            end;
            if(electrode_idx>length(min_idx))
                fprintf('no electrode found!\n');
                flag_found=1;
            end;
        end;
        
        [dummy,min_idx]=sort(electrode_dist_min_com(:));
        for ii=1:3 %show the nearest three contacts
            [ee,cc]=ind2sub(size(electrode_dist_min_com),min_idx(ii));
            fprintf('\tcom dist::<%s_%02d> %2.2f (mm) (%1.1f %1.1f %1.1f)\n',electrode(ee).name,cc,dummy(ii),electrode(ee).coord(cc,1),electrode(ee).coord(cc,2),electrode(ee).coord(cc,3));
        end;
        
        
        
        flag_found=0;
        electrode_idx=1;
        while(~flag_found)
            [ee,cc]=ind2sub(size(electrode_dist_min_com),min_idx(electrode_idx));
            roi(label_idx).electrode_com_dist_electrode_name=electrode(ee).name;
            roi(label_idx).electrode_com_dist_electrode_contact=cc;
            roi(label_idx).electrode_com_dist=dummy(1);
            target_electrode_contact=sprintf('%s%d',electrode(ee).name,cc);
            IndexC = strfind(erf_all(1).name,target_electrode_contact);
            Index = find(not(cellfun('isempty',IndexC)));
            if(~isempty(Index))
                flag_found=1;
                for cond_idx=1:length(erf_all)
                    roi(label_idx).erf_electrode_com_dist(cond_idx).data=squeeze(erf_all(cond_idx).erf_raw(Index,:,:));
                end;
            else
                electrode_idx=electrode_idx+1;
            end;
            
            if(electrode_idx>length(min_idx))
                fprintf('no electrode found!\n');
                flag_found=1;
            end;
        end;
    end;
    
end;

save electrodes_to_labels_061721.mat roi