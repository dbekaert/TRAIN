function [data,format_str_gmt] = load_roipac(datafile,load_flag)
% function which reads roipac data, or data with .rsc file into matlab.
% Currently only supports real 1-band data
% By Bekaert David - Jet Propulsion Laboratory
% modifications:

if nargin <2
    load_flag = 'y';
end


% the rsc file
datafile_rsc = [datafile '.rsc'];
% check if the rsc file exists
if exist(datafile_rsc,'file')~=2
    error([datafile_rsc ' does not exist'])
end

% getting relevant information 
[width] = get_parm_rsc(datafile_rsc,'width');
[length] = get_parm_rsc(datafile_rsc,'file_length');
if isempty(length)
    [length] = get_parm_rsc(datafile_rsc,'length');
end
[type] = get_parm_rsc(datafile_rsc,'format');

%% Getting the precision estimated
fid=fopen(datafile,'r');
fseek(fid, 0, 'eof');
pos=ftell(fid);
byte=pos/width/length;
fseek(fid, 0, 'bof');
switch byte
  case 8 
      format_str_gmt='d';     %double
      format_str='double';     %double
  case 4 
      format_str_gmt='f';     %float, real4
      format_str='float';     %float, real4

  case 2 
      format_str_gmt='i';     %int16
      format_str='int16';     %int16

  case 1 
      format_str_gmt='h';     %int8
      format_str='int8';     %int8
  otherwise
      error('no such format');
end
% give warning as we cannot discriminate between int16 and short
if strcmpi(format_str,'int16')
    fprintf('Could be int16 or short format \n')
end
fclose(fid);


%% in case a precision was given in the file
format_str_given = [];
if ~isempty(type)
    if  strcmpi(type,'r4') || strcmpi(type,'real4')
       format_str_given='float'; 
    elseif strcmpi(type,'h') || strcmpi(type,'short')
        format_str_given='int8';   
        % if the given fortmat is short then update the detection to be short
        if strcmpi(format_str,'int16')
            format_str = 'int8';
        end
    end
else
    fprintf('Will try to figure out precision from the data\n')
end

 
%% checking if the specified precision is different, then that of the automated estimation
if ~isempty(format_str_given)
    if ~strcmpi(format_str_given,format_str)
        fprintf('The automated detected precision is different then what you specified in .rsc file \n:')
        fprintf(['Yours: ' format_str_given '\n'])
        fprintf(['Auto: ' format_str '\n'])

        repeat=1;
        while repeat==1
            action_flag= str2num(input('Keep yours (1), Keep auto (2), Different (3)? [1, 2, or 3] ','s'));
            if isnumeric(action_flag)
                if action_flag==1
                    format_str = format_str_given;
                    repeat=0;
                elseif action_flag==2
                    format_str = format_str;
                    repeat=0;
                elseif action_flag==3
                    action_flag= input('To what do you want to update this? [give a format recognised by GMT]','s');
                    format_str = action_flag;
                    repeat=0;
                end
            end
        end
        fprintf(['Please verify, and re-run'])
    end
end


%% loading the data unless not requested
if strcmpi(load_flag,'y')
    fid = fopen(datafile,'r');
    data = fread(fid,[width, length],format_str);
    fclose(fid);
else
    data = [];
end
