close all; clear all;
clear global  etc_render_fsbrain;

setenv('SUBJECTS_DIR','/Users/fhlin_admin/workspace/seeg/subjects/');
subject='s031';
surf='orig';

file_forward_mat='seeg_fwd_wb_091019.mat';

mri=MRIread('/Users/fhlin_admin/workspace/seeg/subjects/s026/mri/orig.mgz'); %for MAC/Linux

file_stc={
    'seeg_wb_mne_091019_a-vol.stc';
    'seeg_wb_mne_091019_v-vol.stc';
    'seeg_wb_mne_091019_av-vol.stc';
    };

output_stem='';

timeVec=[]; %empty for all time points

target_subject='fsaverage';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load the Talairach transformation matrix from the "pre-OP" data
%talxfm=etc_read_xfm('file_xfm','/Users/fhlin_admin/workspace/seeg/subjects/s026/mri/transforms/talairach.xfm'); %for MAC/Linux

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
            
            fprintf('saving [%s]...\n',fn_output);
            MRIwrite(mri_overlay,fn_output); %native space
            
            if(t_idx==1)
                MRIwrite(mri,fn_under_output); %native space
            end;
            
            %morph
            
            %make a registration file
            fp=fopen('tmp_register.dat','w');
            fprintf(fp,'%s\n',subject);
            fprintf(fp,'1.0\n');
            fprintf(fp,'1.0\n');
            fprintf(fp,'0.15\n');
            fprintf(fp,'1 0 0 0\n');
            fprintf(fp,'0 1 0 0\n');
            fprintf(fp,'0 0 1 0\n');
            fprintf(fp,'0 0 0 1\n');
            fprintf(fp,'round\n');
            fclose(fp);


            if(t_idx==1)
                eval(sprintf('!mri_vol2vol --mov %s --reg tmp_register.dat --o %s%s_under-tal.2mm.mgh  --tal --talres 2',fn_under_output,output_stem,fstem));
                mri_overlay_tal=MRIread(sprintf('%s%s_under-tal.2mm.mgh',output_stem,fstem));
            end;

            %eval(sprintf('!mri_vol2vol --mov %s --reg tmp_register.dat --o %s%s-tal.2mm.mgh  --tal --talres 2',fn_output,output_stem,fstem));
            eval(sprintf('!mri_vol2vol --mov %s --reg tmp_register.dat --o tmp-tal.2mm.mgh  --tal --talres 2',fn_output));

            tmp=MRIread('tmp-tal.2mm.mgh');
            
            mri_overlay_tal.vol(:,:,:,t_idx)=tmp.vol;
        catch ME
        end;
        
        fprintf('\n');
    end;
    
    MRIwrite(mri_overlay_tal,sprintf('%s%s-tal.2mm.mgh',output_stem,fstem));

end;

%without specified electrode
%etc_render_fsbrain('surf','orig','hemi','lh','subject','s026','vol',mri,'talxfm',(talxfm),'overlay_vol_stc',vol_stc,'vol_A',A);

