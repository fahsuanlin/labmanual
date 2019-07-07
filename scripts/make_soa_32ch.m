close all; clear all;

file_soa_mat={
    'VEP_7p5Hz500ms_SMSInIch32_run1-0630-1_sti.mat';
    'VEP_7p5Hz500ms_SMSInIch32_run2-0630-2_sti.mat';
};

output_stem={
    '32ch_soa_01';
    '32ch_soa_02';
    };

for f_idx=1:length(file_soa_mat)
    
    load(file_soa_mat{f_idx});
    
    fstem=output_stem{f_idx};
    file_para=sprintf('%s.para',fstem);
    
    fprintf('writing [%s]....\n',file_para);
    fp=fopen(file_para,'w');
    
    soa_offset=0; %adjustment for slice-timing correction; second
    soa_start=0;
    
    %Cell variable 'onsets' is provided by the matlab file. It includes the SOA
    %(in seconds) of an event
    onset=onsets{1};
    
    for idx=1:length(onset)
        fprintf(fp,'%2.2f\t%d\n',onset(idx)+soa_offset-soa_start,1); %we hard-code the condition '1' to this condition
    end;
    
    fclose(fp);
end;


return;

