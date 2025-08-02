close all; clear all;

t=readtable('Participant_table.csv');

nan_idx=find(isnan(t.MeanTimeToCorrectlyIdentifyMatches_Instance2));

fp=fopen('Participant_table.txt','w');

fprintf(fp,'FID\tIID\tphenotype\n');

for idx=1:size(t,1)
    if(~ismember(idx,nan_idx))
        fprintf(fp,'%d\t%d\t%d\n',t.ParticipantID(idx),t.ParticipantID(idx),t.MeanTimeToCorrectlyIdentifyMatches_Instance2(idx));
    end;
end;
fclose(fp);
