function [] = aps_modis(start_step,end_step)
% [] = aps_modis(start_step,end_step)
% Computes the MODIS delays based on the OSCAR service provied by JPL 
%
%     Copyright (C) 2015  Bekaert David - University of Leeds
%     Email: eedpsb@leeds.ac.uk or davidbekaert.com
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
% By David Bekaert - University of Leeds
% May 2014
%
% Modifications:
% DB	05/2014 	Start from the MERIS scrip the program the function
% DB    07/2014     Have re-visualization of sounding data when already
%                   processed
% DB    08/2014     Include SAR varying conversion factors
% DB    08/2014     Include option to use MODIS recalibrated data
% DB    09/2014     Check if conversion factors are known before
% DB    04/2015     Fix in case filename grd is not starting with OSCAR_Modis*grd
% DB    02/2016     Include aps_systemcall and suppress filename output.



%% loading the StaMPS variables
curdir = pwd;

if start_step==0
   % get the data from the OSCAR webite (JPL) 
   aps_modis_Python
end

if start_step==1 
     % Compute the MERIS conversion factors from the sounding data
    fprintf('\n\nStep 1: Compute the spectrometer conversion factors from the sounding data \n')

    
    % getting parameters from parms_aps file
    sounding_dir = getparm_aps('sounding_dir',1);
    start_date = getparm_aps('sounding_start_date',1);
    end_date = getparm_aps('sounding_end_date',1);    
    time_stamp = getparm_aps('sounding_time_stamp',1);
    sounding_ifg_dates = getparm_aps('sounding_ifg_dates',1);
    start_year_str = start_date(1:4);
    end_year_str = end_date(1:4);
    start_str = start_date(5:6);
    end_str = end_date(5:6);
    time_stamp_str = [];
    for k=1:size(time_stamp,1)
        if k>1
            time_stamp_str = [time_stamp_str '_' time_stamp(k,:)];
        else
            time_stamp_str = [time_stamp(k,:)];
        end
    end
    
     % see if the filename already exist
    if strcmp(sounding_ifg_dates,'y')
        savename = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat'];
    else
        savename = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat'];
    end
    
    if exist(savename,'file')==2
        % check if user wants to visualize or re-run the data
        str = '';
        while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
            str = input(['This file already exist, do you want to re-process the data? [y/n] \n'],'s');
        end

        if strcmpi(str,'y')
            rerun_flag =1;
        else
            sounding_spectrometer_sens_display
            rerun_flag =0;
        end       
    else
       rerun_flag =1;
    end
    
    if rerun_flag ==1;
        sounding_spectrometer_sens;
    end
end

if start_step<=2 && end_step >=2
    % SAR delays using MODIS data
    fprintf('\n\nStep 2: Compute MODIS tropospheric delay for individual dates\n')
    
    % generating a file list of the files to be processed:
    modis_datapath = getparm_aps('modis_datapath',1);
    command_str = 'echo files > modis_batch_file.txt';
    command_str2 = ['ls -d ' modis_datapath filesep '2*' filesep '*.grd >> modis_batch_file.txt'];

    
    aps_systemcall(command_str);
    aps_systemcall(command_str2);
        
    % running the MODIS computation on the file list
    aps_modis_SAR('modis_batch_file.txt')
end
if start_step<=3 && end_step >=3
    fprintf('\n\nStep 3: Computes MODIS tropospheric delay for interferograms\n')
    aps_modis_InSAR
end


cd(curdir)



