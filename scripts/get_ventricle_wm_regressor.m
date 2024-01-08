close all; clear all;

target_subject='s006';
file_register_source={
    '../resting_data/unpack/bold/005/fmcprstc.nii.gz';
    };

file_regression_source={
    '../resting_data/unpack/bold/005/fmcprstc.nii.gz';
    };

file_register={
    './register.dat';
    };

file_aseg={
        'aparc+aseg_fmcprstc_005.nii';
};

file_output={
    'regressor_wm_ventrical_005.mat';
    };
%'regressor_wm_ventricle_123016.mat';


for f_idx=1:length(file_regression_source)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %create a registration matrix from the native subject to the target subject
    %"fsaverage". This registration matrix will be name as
    %"native2fsaverage.dat".
    %eval('!fslregister --s fsaverage --mov bold/004/f.nii --reg ./native2fsaverage.dat --initxfm --maxangle 70');
    if(isempty(file_register{f_idx}))
        [dummy,fstem]=file_output{f_idx};
        file_register{f_idx}=sprintf('%s_register.dat',fstem);
        fprintf('creating registration...[%s]....\r',file_register{f_idx});
        eval(sprintf('!fslregister --s %s --mov %s --reg %s --initxfm --maxangle 70',target_subject, file_register_source{f_idx}, file_register{f_idx}));
    else
        fprintf('using registration...[%s]....\r',file_register{f_idx});
    end;
    
    %apply the "inverse" of the registration such that the aparc+aseg.mgz from
    %"fsaverage" will be transformed to the native subject's anatomical space.
    %The transformed apart+aseg file will be named as "aparc+aseg_f.nii".
    %eval('!mri_vol2vol --mov bold/004/fmc.nii --targ $SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz --reg native2fsaverage.dat --o aparc+aseg_f.nii --inv --interp nearest');
    eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/aparc+aseg.mgz --reg %s --o %s --inv --interp nearest',file_regression_source{f_idx},  target_subject, file_register{f_idx}, file_aseg{f_idx}));
   
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %get white-matter and ventrical from "FreeSurferColorLUT.txt"
    wm=[2, 41];
    % 2   Left-Cerebral-White-Matter              245 245 245 0
    % 41  Right-Cerebral-White-Matter             0   225 0   0
    
    
    ventricle=[4 14 15 43 72 75 76];
    % 4   Left-Lateral-Ventricle                  120 18  134 0
    % 14  3rd-Ventricle                           204 182 142 0
    % 15  4th-Ventricle                           42  204 164 0
    % 43  Right-Lateral-Ventricle                 120 18  134 0
    % 72  5th-Ventricle                           120 190 150 0
    % 75  Left-Lateral-Ventricles                 120 18  134 0
    % 76  Right-Lateral-Ventricles                120 18  134 0
    
    
    d_aseg=MRIread(file_aseg{f_idx});
    v_aseg=d_aseg.vol;
    
    %d_regression=MRIread(file_regression_source);
    %v_regression=d_regression.vol;
    d=MRIread(file_regression_source{f_idx});
    acc=d.vol;
    dim=size(acc);
    %v_regression=reshape(acc,[dim(1)*dim(2)*dim(3), d_regression.nframes]);
    %fprintf('functional data: [%d] voxels x [%d] time points\n',dim(1)*dim(2)*dim(3), d_regression.nframes);
    v_regression=reshape(acc,[dim(1)*dim(2)*dim(3), dim(4)]);
    fprintf('functional data: [%d] voxels x [%d] time points\n',dim(1)*dim(2)*dim(3), dim(4));
    
    tmp=[];
    for wm_idx=1:length(wm)
        idx=find(v_aseg(:)==wm(wm_idx));
        fprintf('wm: index [%d]: [%d] voxels...\n',wm(wm_idx),length(idx));
        tmp=cat(1,tmp,v_regression(idx,:));
    end;
    regressor_wm=mean(tmp,1);
    
    tmp=[];
    for ventricle_idx=1:length(ventricle)
        idx=find(v_aseg(:)==ventricle(ventricle_idx));
        fprintf('ventricle: index [%d]: [%d] voxels...\n',ventricle(ventricle_idx),length(idx));
        tmp=cat(1,tmp,v_regression(idx,:));
    end;
    regressor_ventricle=mean(tmp,1);
    
    save(file_output{f_idx},'regressor_wm','regressor_ventricle');
end;
fprintf('DONE!\n');
