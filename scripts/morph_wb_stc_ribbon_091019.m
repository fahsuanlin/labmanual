close all; clear all;
clear global  etc_render_fsbrain;

setenv('SUBJECTS_DIR','../../subjects/');
subject='s031';
surf='orig';

file_forward_mat='seeg_fwd_wb_091019.mat';

mri=MRIread(sprintf('../../subjects/%s/mri/orig.mgz',subject)); %for MAC/Linux

file_stc={
    'seeg_wb_mne_091019_a_mne-vol.stc';
    'seeg_wb_mne_091019_v_mne-vol.stc';
    'seeg_wb_mne_091019_av_mne-vol.stc';
};

output_stem='';

timeVec=[]; %empty for all time points

target_subject='fsaverage';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the Talairach transformation matrix from the "pre-OP" data
if(ismac)
    talxfm=etc_read_xfm('file_xfm',sprintf('../../subjects/%s/mri/transforms/talairach.xfm',subject)); %for MAC/Linux
    ribbon=MRIread(sprintf('../../subjects/%s/mri/ribbon.mgz',subject));
else
    talxfm=etc_read_xfm('file_xfm',sprintf('../../subjects/%s/mri/transforms/talairach.xfm',subject)); %for MAC/Linux
    ribbon=MRIread(sprintf('../../subjects/%s/mri/ribbon.mgz',subject));
end;

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
    
    if(ismac)
        targ_subj=MRIread('/Applications/freesurfer/average/mni305.cor.subfov2.mgz'); %MNI-Talairach space with 2mm resolution (for MAC)
    elseif(isunix)
        targ_subj=MRIread(sprintf('%s/average/mni305.cor.subfov2.mgz',getenv('FREESURFER_HOME'))); %MNI-Talairach space with 2mm resolution (for server)
    end;
    
    fprintf('loading transformation for subject [%s]...\n',subject);
    targ_xfm=etc_read_xfm('subject',subject);
    
    [overlay_vol, overlay_D,loc_vol_idx]=etc_volstc2vol(vol_stc(:,1:1:end),A,mri,'flag_morph',1,'targ_subj',targ_subj,'targ_xfm',targ_xfm,'vol_ribbon',ribbon);

    fn=sprintf('%s_volstc_tal2mm_ribbon.mgz',fstem);
    fprintf('saving [%s]...\n');
    MRIwrite(overlay_vol,fn); %for MAC/Linux

end;
        
