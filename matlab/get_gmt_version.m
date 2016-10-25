function [gmt5_above, gmt_version]=get_gmt_version()
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


% checking if GMT can be called
% this is a bugfix for the MAC OS systems
command_str = ['gmt --version > gmt_version'];
[gmt_does_not_work, b] = system(command_str);
% Leeds machines define gmt as GMT
if gmt_does_not_work~=0
    command_str = ['GMT --version > gmt_version'];
    [gmt_does_not_work, b] = system(command_str);
end

% GMT does not work well
while gmt_does_not_work~=0
   fprintf('GMT or gmt is not an executable, will try to fix \n')
   % exporting the library path manual
   fprintf('Define the Library path manual in matlab: '''' \n')
   setenv('DYLD_LIBRARY_PATH', '');

   command_str = ['gmt --version > gmt_version'];
   [gmt_does_not_work, b] = system(command_str);
   % Leeds machines define gmt as GMT
   if gmt_does_not_work~=0
       command_str = ['GMT --version > gmt_version'];
       [gmt_does_not_work, b] = system(command_str);
   end

   if gmt_does_not_work==0
       break
   else
       error('Could not fix it, try to fix youself such you can call GMT or gmt from matlab command line as e.g.: gmt --version or GMT --version')
   end
end
clear command_str gmt_does_not_work

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

% trying to call a GMT function 
command_str = ['man psxy > gmt_function_test'];
[gmt_does_not_work, b] = system(command_str);
if gmt_does_not_work~=0
    command_str2 = ['psxy --help> gmt_function_test'];
    [gmt_does_not_work, b] = system(command_str2);
end
%if gmt_does_not_work~=0
command_str2 = ['psxy'];
[gmt_does_not_work, b] = system(command_str2);
if ~isempty(b)
   ix = findstr('usage: psxy',b);
   if ~isempty(ix)
       gmt_does_not_work = 0;
   end
end
%end
% GMT does not work well
while gmt_does_not_work~=0
   fprintf('GMT not an executable, will try to fix \n')
   % exporting the library path manual
   
   fprintf('Define the Library path manual in matlab: '''' \n')
   setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');

   [gmt_does_not_work, b] = system(command_str);
   if gmt_does_not_work==0
       break
   else
       error('Could not fix it, try to fix youself such you can call GMT from matlab command line as e.g.: > gmt --version, and > xyz2grd')
   end
end
clear b command_str gmt_does_not_work
delete('gmt_function_test');


% check if its GMT 5 and above
if str2num(gmt_version(1))>=5
    gmt5_above='y';
else
    gmt5_above='n';
end

% if strcmpi(gmt5_above,'y')
%     fprintf('Will need to use the GMT compatible codes for GMT5 and above... \n')
% else
%     fprintf('Will need to use the GMT compatible codes for version before GMT5... \n')    
% end
