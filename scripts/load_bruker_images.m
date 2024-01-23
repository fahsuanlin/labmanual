function [S,params] = load_bruker_images(varargin)
% S = load_bruker_images(filepath)
%   filepath is the location of the 2dseq file you want to read. 
%
% S = load_bruker_images(filepath,opt1,val1,opt2,val2,...)
%   opt1,val1,opt2,val2 are option name/value pairs.
%
% S = load_bruker_images(opt1,opt1_val,opt2,opt2_val,...)
%   You can also call load_bruker_images without a filepath, in which case
%   a ui file selection box will be presented.
%
% supported options:
%   rescale_recon - true or false (rescales images based on parameters in
%                   the reco file.
%   flip_phase - true or false (are phase & frequency flipped from the
%                standard setup for this image orientation)

if ispc, dirsep = '\';
else dirsep = '/';
end

% linearly rescale images based on recon_slope and offset
rescale_recon = false;
flip_phase = false;
debug_msgs=false;

% if varargin is odd, we assume that the filepath is the first parameter
if(rem(length(varargin),2))
    filepath = varargin{1};
    if(length(varargin)>1)
        varargin = varargin([2:end]);
    else
        varargin = [];
    end
else % otherwise, assume that filepath is empty
    filepath = [];
end

% check options
for i=1:floor(length(varargin)/2)
    option=varargin{i*2-1};
    option_value=varargin{i*2};
    switch lower(option)
    case 'rescale_recon'
    	rescale_recon=option_value;
    case 'flip_phase'
        flip_phase=option_value;     
    case 'debug'
    	debug_msgs=option_value;
    otherwise
        fprintf('unknown option [%s]!\n',option);
        fprintf('error!\n');
        return;
    end;
end;

% parameter values to output
params = [];

% if there was no filepath input, get one using a file selection dialog box
if(nargin==0 || isempty(filepath))
    [filename,imagepath] = uigetfile('*.*', 'Select 2dseq File','2dseq');
    filepath = [imagepath filename];
else
    ind = strfind(filepath,dirsep);
    imagepath = filepath(1:ind(end));
    filename = filepath(ind(end)+1:end);
end

% check the d3proc file
%params.d3proc = read_param_file([imagepath 'd3proc'],debug_msgs);
% lamw: No d3proc file in PV6.0

% check the reco file
params.reco = read_param_file([imagepath dirsep 'reco'],debug_msgs);

% check the method file
params.method = read_param_file([imagepath '..' dirsep '..' dirsep 'method'],debug_msgs);

% check the acqp file
params.acqp = read_param_file([imagepath '..' dirsep '..' dirsep 'acqp'],debug_msgs);

% lamw: Get the matrix size, which used to be in the d3proc file
% check the visu_pars file
params.visu_pars = read_param_file([imagepath dirsep 'visu_pars'],debug_msgs);

% read the image data
if(strcmp(params.reco.RECO_byte_order,'littleEndian'))
    fid = fopen([imagepath dirsep '2dseq'], 'r', 'l');
else
    fid = fopen([imagepath dirsep '2dseq'], 'r', 'b');
end
%if(strcmp(params.method.RECO_wordtype,'_32BIT_SGN_INT'))
% lamw: No RECO_wordtype in method file in PV6.0
if(strcmp(params.reco.RECO_wordtype,'_32BIT_SGN_INT'))
    data = fread(fid,'long');
else
    data = fread(fid,'short');
end
fclose(fid);
if (flip_phase)
    %S = permute(reshape(data,[params.d3proc.IM_SIY params.d3proc.IM_SIX params.d3proc.IM_SIZ params.d3proc.IM_SIT]),[2 1 3 4]);
    % lamw: No d3proc file in PV6.0
    % lamw: Add DIM 5: ACQ_n_movie_frames, which is b value (including b = 0)
    S = permute(reshape(data,[params.visu_pars.VisuCoreSize(2) params.visu_pars.VisuCoreSize(1) params.acqp.NSLICES params.acqp.NR params.acqp.ACQ_n_movie_frames]),[2 1 3 4 5]);
else
    %S = permute(reshape(data,[params.d3proc.IM_SIX params.d3proc.IM_SIY params.d3proc.IM_SIZ params.d3proc.IM_SIT]),[2 1 3 4]);
    % lamw: No d3proc file in PV6.0
    % lamw: Add DIM 5: ACQ_n_movie_frames, which is b value (including b = 0)
%    S = permute(reshape(data,[params.visu_pars.VisuCoreSize(1) params.visu_pars.VisuCoreSize(2) params.visu_pars.VisuCoreFrameCount params.acqp.NSLICES params.acqp.NR params.acqp.ACQ_n_movie_frames]),[2 1 3 4 5 6]);
   % wo: FrameCount = NSLICES? 
%   S = permute(reshape(data,[params.visu_pars.VisuCoreSize(1) params.visu_pars.VisuCoreSize(2) params.acqp.ACQ_n_movie_frames params.acqp.NSLICES]),[2 1 3 4 5 6]);
% wo-use this one   
try
    S = permute(reshape(data,[params.reco.RECO_size(1) params.reco.RECO_size(2) params.reco.RECO_size(3)]),[2 1 3]);
catch
    S = permute(reshape(data,params.reco.RECO_size(1),params.reco.RECO_size(2),[]),[2 1 3]);
end
%    S = permute(reshape(data,[params.visu_pars.VisuCoreSize(1) params.visu_pars.VisuCoreSize(2) params.acqp.NSLICES]),[2 1 3]);
%S = permute(reshape(data,[params.visu_pars.VisuCoreSize(1) params.visu_pars.VisuCoreSize(2)  params.acqp.NSLICES params.acqp.ACQ_n_movie_frames]),[2 1 3 4]);
end

% rescale the reconstruction so that all images have the same window level
if(rescale_recon)
    if(debug_msgs)
        disp('rescale_recon=true: rescaling reconstruction units.');
    end
    %for j=1:params.d3proc.IM_SIZ*params.d3proc.IM_SIT
    % lamw: No d3proc file in PV6.0
%    for j=1:params.acqp.NSLICES*params.acqp.NR*params.acqp.ACQ_n_movie_frames    
% NR & n_movie_frames not workign properly for DTI
    for j=1:params.acqp.NSLICES*params.acqp.ACQ_n_movie_frames
        % WO changed j to 1 in the following line
        S(:,:,j) = (S(:,:,j)-params.reco.RECO_map_offset(1))./params.reco.RECO_map_slope(1);
    end
end

% this function parses the different parameter files
function params = read_param_file(filepath,debug_msgs)
if(exist(filepath,'file')==2)
    if(debug_msgs)
        disp(['reading parameters from ' filepath]);
    end
    fid = fopen(filepath, 'r');
    flag_reading_matrix = false;
    while 1
        tline = fgetl(fid);
        %if strfind(tline,'<5.1>'), disp(filepath), disp(tline), end
        if(flag_reading_matrix)
            if(~ischar(tline) || ~isempty(regexp(tline,'^#')) ...
               || ~isempty(regexp(tline,'^\$')))
                flag_reading_matrix = false;
                
                % if the parameter matrix is numeric, convert it
                %if(~isempty(regexp(matrix,'[\s\d\.]+')) && isnumeric(str2num(matrix)))
                % lamw: isnumeric(str2num(matrix)) thinks <5.1> is numeric,
                %       when it is not, and converts it to NULL. Fix this.
                if(~isempty(regexp(matrix,'[\s\d\.]+')) && ~isempty(str2num(matrix)) && isnumeric(str2num(matrix)))
                    matrix = str2num(matrix);
                end
                if(length(matrix(:))==prod(param_size))
                    % resize according to the param_size, and flip the
                    % first 2 dimensions
                    matrix = permute(reshape(matrix,param_size),[2 1 3:length(param_size)]);
                end
                eval(['params.' param_name '=matrix;']);
            else
                % add elements to the current paramter matrix
                matrix = [matrix char(tline)];
            end
        end
        if ~ischar(tline) break; % end of file
        else
            % if this regexp matches, we are going to be reading in
            % multiple lines into a matrix
            n = regexp(tline,'##\$(.+)=\( (.+) \)$','tokens','once');
            if(~isempty(n))
                flag_reading_matrix = true;
                param_name = char(n(1));
                
                size_str = char(n(2));
                param_size = [];
                while ~isempty(size_str)
                    [t,size_str] = strtok(size_str,', ');
                    param_size = [param_size str2num(char(t))];
                end
                % size must be at least 2 dimension
                if(length(param_size)==1)
                    param_size = [param_size 1];
                end
                % flip dimensions 1 and 2
                param_size = [param_size(2) param_size(1) param_size(3:end)];
                matrix = [];
            else
                % if this regexp matches, we have a single line parameter
                n = regexp(tline,'##\$(.+)=(.+)$','tokens','once');
                if(~isempty(n))
                    param_name = char(n(1));
                    param_value = str2num(char(n(2)));
                    if(isempty(param_value))
                        param_value = char(n(2));
                    end
                    eval(['params.' param_name '=param_value;']); continue;
                end
            end
        end
    end    
    fclose(fid);
else
   error([filepath ' does not exist']);
end
