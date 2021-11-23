close all; clear all;

subject='s031';

path_mri='/Users/fhlin/workspace/seeg/subjects';

file_register='register.dat';

%
xfm=etc_read_xfm('file_xfm',sprintf('%s/%s_post/tmp/%s',path_mri,subject,file_register)); %for MAC/Linux
%xfm=etc_read_xfm('file_xfm','D:\fhlin\Users\fhlin\workspace\seeg\subjects\2036_post\tmp\register.dat'); %for PC

file_mat='electrode_101721_223845.mat'; %electrode coordinates in post-op MRI



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert electrode coordinates from post-op MRI to pre-op MRI

load(file_mat);

electrode_out=electrode;
for e_idx=1:length(electrode)
    
    for c_idx=1:electrode(e_idx).n_contact
        
        surface_coord=electrode(e_idx).coord(c_idx,:);

        surface_coord=inv(xfm)*[surface_coord(:); 1];
        
        electrode_out(e_idx).coord(c_idx,:)=surface_coord(1:3);
        
    end;
end;

electrode=electrode_out;

[dummy,fstem]=fileparts(file_mat);

fprintf('\nmust load [%s] to locate electrodes in ***pre-op*** MRI!\n\n',sprintf('%s_%s.mat',fstem,subject));
save(sprintf('%s_%s.mat',fstem,subject),'electrode');