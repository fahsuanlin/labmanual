close all; clear all;

subject='s031';

path_mri='/Users/fhlin/workspace/seeg/subjects';

file_register='register.dat';

targ=MRIread('/Applications/freesurfer/subjects/fsaverage/mri/orig.mgz');

file_mat='electrode_101721_223845.mat'; %electrode coordinates in post-op MRI

file_roi='electrodes_to_labels_061721.mat';
roi_index=[1 4]; %auditory cortex; left and right hemisphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert electrode coordinates from post-op MRI to MNI MRI
mri=MRIread(sprintf('%s/%s/mri/orig.mgz',path_mri,subject)); %for MAC/Linux
%
xfm=etc_read_xfm('file_xfm',sprintf('%s/%s_post/tmp/%s',path_mri,subject,file_register)); %for MAC/Linux
%xfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036_post\tmp\register.dat'); %for PC

fprintf('loading transformation for subject %s]...\n',subject);
mov_xfm=etc_read_xfm('subject',subject);

load(file_mat);

electrode_out=electrode;
for e_idx=1:length(electrode)
    
    for c_idx=1:electrode(e_idx).n_contact
        
        surface_coord=electrode(e_idx).coord(c_idx,:);

        surface_coord=targ.tkrvox2ras*inv(targ.vox2ras)*mov_xfm*mri.vox2ras*inv(mri.tkrvox2ras)*inv(xfm)*[surface_coord(:); 1];
        
        electrode_out(e_idx).coord(c_idx,:)=surface_coord(1:3);
        
    end;
end;

electrode=electrode_out;

[dummy,fstem]=fileparts(file_mat);

fprintf('\nmust load [%s] to locate electrodes in MNI MRI!\n\n',sprintf('%s_tal_%s.mat',fstem,subject));

save(sprintf('%s_tal_%s.mat',fstem,subject),'electrode');



for roi_idx=1:length(roi_index)
    
    %find electrode and contact index
    if(~isempty(file_roi))
        load(file_roi);
        if(~isempty(electrode))
            for electrode_idx=1:length(electrode)
                if(strcmp(electrode(electrode_idx).name,roi(roi_index(roi_idx)).electrode_min_dist_electrode_name))
                    electrode_idx_select=electrode_idx;
                end;
            end;
        end;        
        contact_idx_select=roi(roi_index(roi_idx)).electrode_min_dist_electrode_contact;
    else
        roi=[];
    end;
    
    E=electrode(electrode_idx_select);
    E.name=sprintf('%s_%s',E.name,subject);
    E.contact_idx=contact_idx_select;
    
    save(sprintf('electrode_tal_mri_%s_061721.mat',roi(roi_index(roi_idx)).name),'E');
    
end;
