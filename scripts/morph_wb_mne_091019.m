close all; clear all;
clear global  etc_render_fsbrain;

%setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/seeg/subjects/');
setenv('SUBJECTS_DIR','/space/maki7/1/fhlin/seeg/subjects/');
subject='s031';
surf='orig';

file_forward_mat='seeg_fwd_wb_091019.mat';

%mri=MRIread('/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/orig.mgz'); %for MAC/Linux
mri=MRIread('/space/maki7/1/fhlin/seeg/subjects/s031/mri/orig.mgz'); %for MAC/Linux

file_stc={
%    'seeg_wb_mne_091019_a_mne-vol.stc';
%    'seeg_wb_mne_091019_v_mne-vol.stc';
%    'seeg_wb_mne_091019_av_mne-vol.stc';
    'seeg_wb_mne_091019_a-vol.stc';
    'seeg_wb_mne_091019_v-vol.stc';
    'seeg_wb_mne_091019_av-vol.stc';
    };

output_stem='';

timeVec=[]; %empty for all time points

target_subject='fsaverage';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load the Talairach transformation matrix from the "pre-OP" data
%talxfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/seeg/subjects/s031/mri/transforms/talairach.xfm'); %for MAC/Linux

load(file_forward_mat);
%append the forward matrix object 'A' surface coordinates; used for
%interpolation
hemi={'lh','rh'};
for hemi_idx=1:2
    file_surf=sprintf('%s/%s/surf/%s.%s',getenv('SUBJECTS_DIR'),subject,hemi{hemi_idx},surf);
    [A(hemi_idx).vertex_coords, A(hemi_idx).faces] = read_surf(file_surf);
end;


for f_idx=1:length(file_stc)
    
    [dummy, fstem]=fileparts(file_stc{f_idx});
    
    fn_output=sprintf('%s%s.mgz',output_stem,fstem);
    fn_under_output=sprintf('%s%s_under.mgz',output_stem,fstem);
    
    %volume stc
    [vol_stc,d0,d1]=inverse_read_stc(file_stc{f_idx});
    
    overlay_D{1}=[];
    overlay_D{2}=[];
    loc_vol_idx{1}=[];
    loc_vol_idx{2}=[];
    
    
    mri_overlay=mri;
    
    if(isempty(timeVec))
        timeVec=[1:size(vol_stc,2)];
    end;
    
    for t_idx=1:length(timeVec)
        
        fprintf('time point [%04d]...',timeVec(t_idx));
        
        try
            %initialize
            loc_vol=[];
            for hemi_idx=1:2
                
                n_dip(hemi_idx)=size(A(hemi_idx).A,2);
                n_source(hemi_idx)=n_dip(hemi_idx)/3;
                
                switch hemi_idx
                    case 1
                        offset=0;
                        hemi_str='lh';
                    case 2
                        offset=n_source(1);
                        hemi_str='rh';
                end;
                
                %get source estimates at cortical and sub-cortical locations
                X_hemi_cort=vol_stc(offset+1:offset+length(A(hemi_idx).v_idx),t_idx);
                X_hemi_subcort=vol_stc(offset+length(A(hemi_idx).v_idx)+1:offset+n_source(hemi_idx),t_idx);
                
                
                %smooth source estimates at cortical locations
                ov=zeros(size(A(hemi_idx).vertex_coords,1),1);
                ov(A(hemi_idx).v_idx+1)=X_hemi_cort;
                
                flag_overlay_D_init=1;
                %if(isfield(etc_render_fsbrain,'overlay_D'))
                %    if(length(etc_render_fsbrain.overlay_D)==2)
                %        if(~isempty(etc_render_fsbrain.overlay_D{hemi_idx}))
                if(~isempty(overlay_D{hemi_idx}))
                    flag_overlay_D_init=0;
                end;
                %        end;
                %    end;
                %end;
                %if(flag_overlay_D_init) etc_render_fsbrain.overlay_D{hemi_idx}=[];end;
                
                overlay_smooth=5;
                fprintf('smoothing <%s>...',hemi_str);
                [ovs,dd0,dd1,overlay_Ds,overlay_D{hemi_idx}]=inverse_smooth('','vertex',A(hemi_idx).vertex_coords','face',A(hemi_idx).faces','value',ov,'value_idx',A(hemi_idx).v_idx+1,'step',overlay_smooth,'n_ratio',length(ov)/size(X_hemi_cort,1),'flag_display',0,'flag_regrid',0,'flag_fixval',0,'D',overlay_D{hemi_idx});
                
                %assemble smoothed source at cortical locations and sources at sub-cortical locations
                X_wb{hemi_idx}=cat(1,ovs(:),X_hemi_subcort(:));
                
                flag_cal_loc_vol_idx=1;
                %if(isfield(etc_render_fsbrain,'loc_vol_idx'))
                %    if(length(etc_render_fsbrain.loc_vol_idx)==2)
                %        if(~isempty(etc_render_fsbrain.loc_vol_idx{hemi_idx}))
                if(~isempty(loc_vol_idx{hemi_idx}))
                    flag_cal_loc_vol_idx=0;
                end;
                %        end;
                %    end;
                %end;
                
                fprintf('creating indices <%s>...',hemi_str);
                
                if(flag_cal_loc_vol_idx==1)
                    etc_render_fsbrain.loc_vol_idx{hemi_idx}=[];
                    %get coordinates from surface to volume
                    loc=cat(1,A(hemi_idx).vertex_coords./1e3,A(hemi_idx).wb_loc);
                    loc_surf=[loc.*1e3 ones(size(loc,1),1)]';
                    tmp=inv(mri.tkrvox2ras)*loc_surf;
                    loc_vol{hemi_idx}=round(tmp(1:3,:))';
                    loc_vol{hemi_idx}=round(tmp(1:3,:))';
                    loc_vol_idx{hemi_idx}=sub2ind(size(mri.vol),loc_vol{hemi_idx}(:,2),loc_vol{hemi_idx}(:,1),loc_vol{hemi_idx}(:,3));
                end;
                
                
            end;
            
            
            fprintf('saving data ...');
            tmp=zeros(size(mri.vol));
            
            for hemi_idx=1:2
                tmp(loc_vol_idx{hemi_idx})=X_wb{hemi_idx};
            end;
            
            mri_overlay.vol=tmp;
        
            if(t_idx==1)
                if(strcmp(target_subject,'fsaverage'))
                    %targ=MRIread('/Applications/freesurfer/average/mni305.cor.subfov2.mgz'); %MNI-Talairach space with 2mm resolution (for MAC)
                    targ=MRIread(sprintf('%s/average/mni305.cor.subfov2.mgz',getenv('FREESURFER_HOME'))); %MNI-Talairach space with 2mm resolution (for server)

                    fprintf('loading transformation for subject p%s]...\n',subject);
                    mov_xfm=etc_read_xfm('subject',subject);
                end;
                
                
                %mov=MRIread(fn_under_output);
                R=inv(mri.vox2ras)*inv(mov_xfm)*(targ.vox2ras);
                MRI_underlay_tal=etc_MRIvol2vol(mri,targ,R);
                MRIwrite(MRI_underlay_tal,sprintf('%s%s_under-tal.2mm.mgh',output_stem,fstem));
                
                mri_overlay_tal=MRI_underlay_tal;
                mri_overlay_tal.nframes=length(timeVec);
                mri_overlay_tal.vol=zeros(mri_overlay_tal.volsize(1),mri_overlay_tal.volsize(2),mri_overlay_tal.volsize(3),length(timeVec));
            end;
	    %mov=MRIread(fn_output);
            R=inv(mri_overlay.vox2ras)*inv(mov_xfm)*(targ.vox2ras);
            tmp=etc_MRIvol2vol(mri_overlay,targ,R);
            mri_overlay_tal.vol(:,:,:,t_idx)=tmp.vol;
        catch ME
        end;
        
        fprintf('\n');
    end;
    
    MRIwrite(mri_overlay_tal,sprintf('%s%s-tal.2mm.mgh',output_stem,fstem));

end;

%without specified electrode
%etc_render_fsbrain('surf','orig','hemi','lh','subject','s031','vol',mri,'talxfm',(talxfm),'overlay_vol_stc',vol_stc,'vol_A',A);

