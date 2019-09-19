close all; clear all;

file_soa='SOA_s026.mat';

soa_var{1}={'R1_Onset_A','R1_Onset_V','R1_Onset_VA','R1_Onset_AM','R1_Onset_VM', 'R1_Onset_VAM'};
soa_code{1}=[1 2 3 11 21 31];
soa_var{2}={'R2_Onset_A','R2_Onset_V','R2_Onset_VA','R2_Onset_AM','R2_Onset_VM', 'R2_Onset_VAM'};
soa_code{2}=[1 2 3 11 21 31];

output_stem='fmri_soa';

load(file_soa);

for f_idx=1:length(soa_var)
    
    %load(file_soa_mat{f_idx});
    for ii=1:length(soa_var{f_idx})
        eval(sprintf('tmp.time=%s;',soa_var{f_idx}{ii}));
        tmp.event=ones(size(tmp.time)).*soa_code{f_idx}(ii);
        if(ii==1)
            trigger=tmp;
        else
            trigger=etc_trigger_append(trigger,tmp);
        end;
    end;
    
    fstem=output_stem;
    file_para=sprintf('%s_%02d.para',fstem,f_idx);
    
    fprintf('writing [%s]....\n',file_para);
    fp=fopen(file_para,'w');
    
    for idx=1:length(trigger.time)
        fprintf(fp,'%2.2f\t%d\n',trigger.time(idx),trigger.event(idx));
    end;
    
    fclose(fp);
end;