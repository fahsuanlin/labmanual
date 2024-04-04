close all; clear all;

%addpath('/autofs/space/maki6_001/users/eva/toolbox/tool_mb_recon/VD/');
%addpath('/autofs/space/maki6_001/users/eva/toolbox/tool_mb_recon/');

addpath('/home/fhlin/matlab/toolbox/fhlin_toolbox/codes/tool_mb_recon/VD');
addpath('/home/fhlin/matlab/toolbox/fhlin_toolbox/codes/tool_mb_recon/');

path_dat_source={
%'/space/maki6/1/users/fhlin/eegmri/180330_SYH/fmri_raw/meas/';
'/space_lin2/fhlin/eegmri/180330_SYH/fmri_raw/meas/';
};
path_mat_destination='/space/maki6/1/users/fhlin/eegmri/180330_SYH/fmri_analysis';


for f_idx=1:length(path_dat_source)
        ref_files=dir(sprintf('%s/*MBSIREPI*ref*.dat',path_dat_source{f_idx}));
        acc_files=dir(sprintf('%s/*MBSIREPI*acc*.dat',path_dat_source{f_idx}));

        for i=1:length(acc_files)

                [ ref,EPInew ] = MBrecon( [path_dat_source{f_idx},ref_files(i).name],[path_dat_source{f_idx},acc_files(i).name],0.0005 );
                %acc=abs(EPInew(:,:,:,61:end));
                acc=abs(EPInew);

                save(sprintf('%s/mb_run_%s_ref.mat',path_mat_destination,num2str(i)),'ref','-v7.3');
                save(sprintf('%s/mb_run_%s_acc.mat',path_mat_destination,num2str(i)),'acc','-v7.3');
        end
end;
