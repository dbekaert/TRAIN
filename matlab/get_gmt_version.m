function [gmt5_above, gmt_version,GMT_string]=get_gmt_version()
% function that checks the GMT version and returns version number and if
% its GMT5 and above. This function also includes a fix for MAC OS. 
%
%     Copyright (C) 2016  Bekaert David 
%     davidbekaert.com
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
% By Bekaert David
% modifications
% 02/2016   DB      Fix in case the dump of GMT --version does not work
% 03/2016   DB      Fix in case gmt is defined as GMT executable
% 03/2016   DB      Be more clear on the GMT and gmt definition
% 07/2016   DB      Try to catch the case when the GMT manual was not installed correctly
% 08/2016   DB      Try cathing in case the library path was not set properly
% 06/2017   DB      Add linux library path, add gmt string to be preceded by gmt functions in matlab.


% checking if gmt can be called and whether it is gmt or GMT
TRY_gmt_executables{1}='gmt ';
TRY_gmt_executables{2}='GMT ';
counter =1;
gmt_does_not_work=1;
GMT_string=[];
while gmt_does_not_work~=0 && counter<=length(TRY_gmt_executables)
    command_str = [TRY_gmt_executables{counter} ' --version > gmt_version'];
    [gmt_does_not_work, b] = system(command_str);
    
    % break in case it works
    if gmt_does_not_work==0
        fprintf(['GMT works with executable: ' TRY_gmt_executables{counter} '\n' ]);
        GMT_string = TRY_gmt_executables{counter};
        break
    end
    counter = counter+1;
end
clear counter


% if it did not work try to see if the library path can be changed to find GMT instead.
if gmt_does_not_work~=0 
   fprintf('GMT or gmt is not an executable, will try to fix the library paths \n')
   counter2 =1;
   counting =1;
   TRY_PATHS{1} = 'DYLD_LIBRARY_PATH';      % typical for mac
   TRY_PATHS{2} = 'LD_LIBRARY_PATH';        % typical for linux
   TRY_SET_PATHS{1}='';                     % try resetting
   TRY_SET_PATHS{2}='/usr/local/bin/';      % try defining usr local bin
   TRY_SET_PATHS{3}='/usr/bin/';            % try defining usr bin
   
   while gmt_does_not_work~=0 && counting<=length(TRY_PATHS)*length(TRY_SET_PATHS)
       % exporting the library path manual
       for counter3=1:length(TRY_SET_PATHS)
           fprintf(['Reset the Library path in matlab: ' TRY_PATHS{counter2} '=''' TRY_SET_PATHS{counter3}    ''' \n'])
           setenv(TRY_PATHS{counter2}, TRY_SET_PATHS{counter3});

           % loop over the GMT or gmt option as executable
           counter=1;
           while gmt_does_not_work~=0 && counter<=length(TRY_gmt_executables)
                command_str = [TRY_gmt_executables{counter} ' --version > gmt_version'];
                [gmt_does_not_work, b] = system(command_str);
                % break in case it works
                if gmt_does_not_work==0
                    fprintf(['GMT works with executable: ' TRY_gmt_executables{counter} '\n' ]);
                    GMT_string = TRY_gmt_executables{counter};
                    break
                end
                counter = counter+1;
           end
           clear counter command_str 
       
      
           % report back
           if gmt_does_not_work==0
               break
           elseif gmt_does_not_work==1 & counter2==length(TRY_PATHS) & counter3==length(TRY_SET_PATHS)
               error('Could not fix it, try to fix youself such you can call GMT or gmt from matlab command line as e.g.: gmt --version or GMT --version')
           end
           counter3 = counter3+1;
           counting = counting+1;

       end
       counter2 = counter2+1;
    end
end


% loading the version number
fid = fopen('gmt_version');
gmt_version = textscan(fid,'%s');
gmt_version = gmt_version{1};
fclose(fid);

% somtimes the machine does not want to dump the information lets try to
% retrieve different
if isempty(gmt_version)
   if ~isempty(b)
       ix = findstr('GMT Version',b);
       if ~isempty(ix)
           gmt_version =strtrim(b(ix+12-1:ix+12-1+2));
       else
          error('Could not retrieve your GMT version') 
       end
   end
else
    gmt_version =gmt_version{:};
    fprintf(['You are using GMT version: ' gmt_version '\n']);
end
clear b
delete('gmt_version');


% check if the gmt function is an executatble on its own or needs to have GMT/gmt in front of it
gmt_function_does_not_work=1;
command_str = ['psxy'];
[gmt_function_does_not_work, b] = system(command_str);
if ~isempty(b)
   ix = findstr('usage: psxy',b);
   if ~isempty(ix)
       gmt_function_does_not_work = 0;
       % no gmt or GMT is needed in front of the command
       GMT_string = [];
   end
end
% if not check if the function can be found by putting gmt in fromt
if gmt_function_does_not_work~=0
    [gmt_function_does_not_work, b] = system([GMT_string ' ' command_str]);
    if ~isempty(b)
       ix = findstr('usage: psxy',b);
       if ~isempty(ix)
           gmt_function_does_not_work = 0;
           fprintf(['WARNING: GMT functions need to be preceded by ' GMT_string(1:end-1) ': e.g. ' GMT_string 'xyz2grd \n'])
       end
    end
end
% GMT does not work well
if gmt_function_does_not_work~=0
    error('Could not find the gmt functions as executable (e.g. xyz2grd) or as gmt function (gmt xyz2grd) executable. Try yourself on commandline ')
end


% check if its GMT 5 and above
if str2num(gmt_version(1))>=5
    gmt5_above='y';
else
    gmt5_above='n';
end
