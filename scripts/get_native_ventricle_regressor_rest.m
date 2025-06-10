close all; clear all;

root_path='/space_lin1/hcp';
root_path='/Users/fhlin/workspace/hcp'

target_subject='fsaverage';

file_aseg='aparc+aseg_f.nii';

setenv('SUBJECTS_DIR','/space_lin1/hcp/subjects');
setenv('SUBJECTS_DIR','/Users/fhlin/workspace/hcp/subjects');

fstem={
'rfMRI_REST1_LR_hp2000_clean';
%'rfMRI_REST1_RL_hp2000_clean';
%'rfMRI_REST2_LR_hp2000_clean';
%'rfMRI_REST2_RL_hp2000_clean';
};

output_fstem={
'native_regressors_rest1_LR';
%'native_regressors_rest1_RL';
%'native_regressors_rest2_LR';
%'native_regressors_rest2_RL';
}

error_subject={};

subject='';
d=textread('subject_list_all.txt');
%for d_idx=1:length(d)
for d_idx=1:1
        subject{d_idx}=num2str(d(d_idx));
end;



for stem_idx=1:length(fstem)
for f_idx=1:length(subject)
	fprintf('[%s]...[%04d||%04d]...\r',subject{f_idx},f_idx,length(subject));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %create a registration matrix from the native subject to the target subject
    %"fsaverage". This registration matrix will be name as
    %"native2fsaverage.dat".
    %eval('!fslregister --s fsaverage --mov bold/004/f.nii --reg ./native2fsaverage.dat --initxfm --maxangle 70');
    %eval(sprintf('!fslregister --s %s --mov %s --reg %s --initxfm --maxangle 70',target_subject, file_register_source{f_idx}, file_register));
    
    %apply the "inverse" of the registration such that the aparc+aseg.mgz from
    %"fsaverage" will be transformed to the native subject's anatomical space.
    %The transformed apart+aseg file will be named as "aparc+aseg_f.nii".
    %eval('!mri_vol2vol --mov bold/004/fmc.nii --targ $SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz --reg native2fsaverage.dat --o aparc+aseg_f.nii --inv --interp nearest');

    %registration has been done previously
    %file_register=sprintf('/space_lin1/hcp/analysis/%s_register.dat',subject{f_idx});
    file_register=sprintf('/Users/fhlin/workspace/hcp/analysis/%s_register.dat',subject{f_idx});
    
    file_mov=sprintf('%s/%s/Preprocessed/%s.nii.gz',root_path,subject{f_idx},fstem{stem_idx});
   
    if(exist(file_mov))

    file_output=sprintf('/space_lin1/hcp/analysis/%s_%s.mat',subject{f_idx},output_fstem{stem_idx});
    if(~exist(file_output))
%    eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/aparc+aseg.mgz --reg %s --o %s --inv --interp nearest',file_mov,  target_subject, file_register, file_aseg));
   
	try 
	 %eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/aparc+aseg.mgz --reg %s --o /space_lin1/hcp/analysis/%s_%s --inv --interp nearest',file_mov,  subject{f_idx}, file_register, subject{f_idx}, file_aseg));
   eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/aparc+aseg.mgz --reg %s --o /Users/fhlin/workspace/hcp/analysis/%s_%s --inv --interp nearest',file_mov,  subject{f_idx}, file_register, subject{f_idx}, file_aseg));
   
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
    
    
    %d_aseg=MRIread(sprintf('/space_lin1/hcp/analysis/%s_%s',subject{f_idx},file_aseg));
    d_aseg=MRIread(sprintf('/Users/fhlin/workspace/hcp/analysis/%s_%s',subject{f_idx},file_aseg));
    v_aseg=d_aseg.vol;
    
    %d_regression=MRIread(file_regression_source);
    %v_regression=d_regression.vol;
    d=MRIread(file_mov);
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
    
    %file_output=sprintf('/space_lin1/hcp/analysis/%s_%s.mat',subject{f_idx},output_fstem{stem_idx});
    file_output=sprintf('/Users/fhlin/workspace/hcp/analysis/%s_%s.mat',subject{f_idx},output_fstem{stem_idx});

    save(file_output,'regressor_wm','regressor_ventricle');
    
    %eval(sprintf('!rm /space_lin1/hcp/analysis/%s_%s',subject{f_idx},file_aseg'));


	catch
		error_subject{end+1}=subject{f_idx};
	end;
    else
	fprintf('regressor [%s] existed!\n',file_output);
    end;
    end;
end;
fprintf('\n');

fprintf('error subject:\n');
for f_idx=1:length(error_subject)
fprintf('[%s]\n',error_subject{f_idx});
end;
end;

fprintf('DONE!\n');
