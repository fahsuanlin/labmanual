% Extract Bruker scan details and write to file

%clear

if ispc, dirsep = '\';
else     dirsep = '/';
end

%data_dir = 'C:\Users\woakden\Desktop\CEST\JB1-Week4\20161020_084058_09_16_JB1_C0M0_1_3';
data_dir = pwd;

if ~exist(data_dir,'dir'), error('Directory not found.'), end
    
%% Read comment and time_run in procpar files

D = dir(data_dir);
sercount = 0; % Series count
for i = 1:numel(D)
    if D(i).isdir && exist([data_dir dirsep D(i).name dirsep 'pdata' dirsep '1' dirsep '2dseq'],'file')
        sercount = sercount + 1;
        expnum(sercount) = str2double(D(i).name);
    end
end
expnum = sort(expnum);

scan_details = cell([sercount 16]);
for i = 1:sercount
    disp(expnum(i))
    debug_msgs = false;
    % Read method file
    method = read_param_file([data_dir dirsep num2str(expnum(i)) dirsep 'method'],debug_msgs);
    % Read acqp file
    acqp = read_param_file([data_dir dirsep num2str(expnum(i)) dirsep 'acqp'],debug_msgs);
    % Read reco file
    reco = read_param_file([data_dir dirsep num2str(expnum(i)) dirsep 'pdata' dirsep '1' dirsep 'reco'],debug_msgs);
    


    try
        scan_details{i, 1} = strrep(datestr(strrep(acqp.ACQ_time(2:20),'T',' '),30),'T',' ');
    catch ME
        scan_details{i, 1} = '0';
    end
        scan_details{i, 2} = expnum(i);
        scan_details{i, 3} = acqp.ACQ_scan_name(2:end-1);
        scan_details{i, 4} = method.PVM_RepetitionTime;
        if ((strcmp(method.Method,'<Bruker:POSITION>'))), scan_details{i, 5} = 0;
        else
            if (isfield(method, 'PVM_EchoTime')), scan_details{i, 5} = method.PVM_EchoTime;
            else  scan_details{i, 5} = method.EchoTime; 
            end
        end

        scan_details{i, 6} = acqp.ACQ_flip_angle;
        
        if ((strcmp(method.Method,'<Bruker:STEAM>'))||(strcmp(method.Method,'<Bruker:PRESS>'))||(strcmp(method.Method,'<User:mtPRESS>'))), 
            scan_details{i,7} = method.PVM_VoxArrSize(1); 
            scan_details{i,8} = method.PVM_VoxArrSize(2); 
            scan_details{i,9} = method.PVM_VoxArrSize(3); 
            scan_details{i,10} = 1;
            scan_details{i,11} = 1;

        else
            scan_details{i,7} = method.PVM_Fov(1);
           try
                scan_details{i, 8} = method.PVM_Fov(2);
            catch ME
                scan_details{i, 8} = 0;
            end
            scan_details{i, 9} = method.PVM_SliceThick;
            scan_details{i,10} = method.PVM_Matrix(1);
            try
                scan_details{i,11} = method.PVM_Matrix(2);
            catch ME
                scan_details{i, 8} = 0;
            end
            try
                scan_details{i,12} = method.PVM_Matrix(3); % 3D scan
            catch ME
                scan_details{i,12} = acqp.NSLICES;  % 2D scan
            end
            %if strfind(method.PVM_SpatDimEnum,'2D'), scan_details{i,12} = acqp.NSLICES;  % 2D scan
            %else scan_details{i,12} = method.PVM_Matrix(3); % 3D scan
            %end
        end
        
        scan_details{i,13} = method.PVM_DigSw;
        
        try 
            scan_details{i,14} = method.PVM_NRepetitions;
        catch ME
            scan_details{i,14} = 1;
        end
        %if strcmp(method.Method,'<Bruker:FieldMap>'), scan_details{i,3} = '<Bruker:FieldMap>'; scan_details{i,14} = 1;
        %else scan_details{i,14} = method.PVM_NRepetitions;
        %end
        scan_details{i,15} = acqp.RG;
        if (isfield(method, 'PVM_SelIrInvTime')), scan_details{i,16} = method.PVM_SelIrInvTime;
        else scan_details{i,16} = 0; end

        if isfield(method,'PVM_MapShimStatus'), scan_details{i,17} = method.PVM_MapShimStatus; end % Not in PV5.1
        %if strcmp(method.Method,'<Bruker:FieldMap>'), scan_details{i,18} = 1;
        %else
            scan_details{i,18} = method.PVM_NAverages;
    try
        scan_details{i, 19} = acqp.ACQ_time(13:20);
    catch ME
        scan_details{i, 19} = '0';
    end
    try
        scan_details{i, 20} = reco.RECO_time(13:20);
    catch ME
        scan_details{i, 20} = '0';
    end

end
scan_details = sortrows(scan_details,1); % Sort by date and time

%% Write to file

fid = fopen([data_dir dirsep 'scan details.txt'],'w');
fprintf(fid,'Exp        Series Name    TR/ TE     FA   FOV(Re*Ph*Sl)   Matrix  Slices   BW(Re)  Reps  RecGain TI Start Time End Time\n');
for i = 1:sercount
    if scan_details{i, 8} < 10
        fprintf(fid,'%2.f %18s  TR/TE:%4.f/%3.fms FA:%3.f° FOV:%3.fx%3.fx%2.1fmm  %3.fx%3.fx%3.f BW:%6.fHz NEX:%1.f RG:%3.f TI:%4.fms  %s %s\n',scan_details{i,2:16},scan_details{i,19:20});
    else
        fprintf(fid,'%2.f %18s  TR/TE:%4.f/%3.fms FA:%3.f° FOV:%3.fx%3.fx%2.f mm  %3.fx%3.fx%3.f BW:%6.fHz NEX:%1.f RG:%3.f TI:%4.fms  %s %s\n',scan_details{i,2:16},scan_details{i,19:20});
    end
end
fclose(fid);

%% Write to file 2

fid = fopen([data_dir dirsep 'scan timing.txt'],'w');
fprintf(fid,'Exp        Series Name    Start Time End Time\n');
for i = 1:sercount
    if scan_details{i, 8} < 10
        fprintf(fid,'%2.f %18s  %s %s\n',scan_details{i,2:3},scan_details{i,19:20});
    else
        fprintf(fid,'%2.f %18s  %s %s\n',scan_details{i,2:3},scan_details{i,19:20});
    end
end
fclose(fid);