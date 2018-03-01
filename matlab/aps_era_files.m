function [] = aps_era_files(orderflag_ECMWF_website,era_model)
% script that runs and checks which ERA-I or ERA5 data files are needed based on
% the satellite pass time from the ECMWF website.
% INPUTS:
%
%        ****** Specific for ECMWF data website ****** NOT BADC
%        orderflag_ECMWF_website = 0/1 - no request(0)  (default(0)),
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
% By David Bekaert - October 2013
% University of Leeds
% Modifications:
% HB 01/2014    adding automatic reading of UTC_sat & ordering of necessary
%               era files from ECMWF data website
% DB 02/2014    Include non-stamps processed data support 
% DB 02/2014    Put the python path in the source_file
% DB 03/2014    Change script such it auto downloads the data from ECMWF website
% DB 03/2014    Suppress command window output
% DB 06/2014    Include an overwrite flag for files and delete the files
%               first otherwize download does not work
% DB 06/2014    Fixed typo in the East extend of the ERA data downlaod
% DB 07/2014    Redefine meris_lat(lon)_range to region_lat(lon)_range
% DB 06/2015    Fix typo
% DB 11/2017    Stagger downloads by 5 sec as new api has a maximum download.
% DB 02/2018    Allow for acquisition variable UTC times, integrate ERA 5


stdargin = nargin ; 
if stdargin<1
    orderflag_ECMWF_website = 0;
end
if nargin<2
    era_model='era';
end

% getting the variables from the parms_aps file
workdir = pwd;

% datapath specific
era_datapath = getparm_aps('era_datapath');
if isempty(era_datapath)
    error('please specify era_datapath')
end

% loading the data
stamps_processed = getparm_aps('stamps_processed');
UTC_sat =  getparm_aps('UTC_sat');
ifgday_matfile = getparm_aps('ifgday_matfile');
ifgs_dates = load(ifgday_matfile);
ifg_dates = ifgs_dates.ifgday;
if ~strcmp(stamps_processed,'y')
    ifg_dates = datenum(num2str(ifg_dates),'yyyymmdd');
end


%% Compute based on satellite pass which weather model outputs that will be used
[time_before,time_after, date_before, date_after,f_before,f_after] = aps_weather_model_times(era_model,ifg_dates,UTC_sat);
time_vector = [time_before ; time_after];
date_vector = [date_before ;  date_after];
clear f_before f_after


%% generate a loop over the number of unique SAR acquistions
filelist = [];
datelist = [];
n_dates = size(date_vector,1);

filelist = [repmat('ggap',n_dates,1) date_vector time_vector repmat('.nc',n_dates,1)];         %Format ggapYYYYMMDDHHMM.nc;
datelist = date_vector;                                                                                  %Format YYYYMMDD;
filelist = unique(filelist,'rows');
datelist = unique(datelist,'rows');

% outputing this information to a file
if strcmpi(era_model,'era')
    fid = fopen('ERA_I_files.txt','w');
    fid2 = fopen('ERA_I_dates.txt','w');
elseif strcmpi(era_model,'era5')
    fid = fopen('ERA5_files.txt','w');
    fid2 = fopen('ERA5_dates.txt','w');
end 
fprintf(['Required ' era_model ' files for all interferograms \n'])
for k=1:size(filelist,1)
    fprintf([filelist(k,:) '\n']);
    fprintf(fid,[filelist(k,:) '\n']);
end
fclose(fid);

for k=1:size(datelist,1)
    fprintf(fid2,[datelist(k,3:end) '\n']);
end
fclose(fid2);



%% Below is specific for the ECMWF data website
crop_range_in = 2; % increasing extent of weather data region by this value (degree)
overwrite_flag=-1;
if orderflag_ECMWF_website==1
    % extend the BBOX a bit larger than the data bbox
    [S,N,W,E] = aps_weather_model_crop;
    
    % weather model region
    fprintf('getting the data from the ECMWF service ...')
    weatherregstr = [N,'/',W,'/',S,'/',E];   % N/W/S/E
    fprintf('weather model region N/W/S/E %s \n',weatherregstr);
    fprintf('using mars service from ECMWF downloading to \n %s \n',era_datapath);
    
    for l = 1:size(filelist,1)
        subdirpath = [era_datapath,filesep,filelist(l,5:12),filesep];
        % make date dir if needed
        if exist(subdirpath,'dir')~=7
            mkdir(subdirpath)
        end
        % check if the ECMWF file exists
        if exist([subdirpath,filelist(l,:)],'file') == 0
            cd(subdirpath)
            fprintf('Order and downloading %s \n',filelist(l,5:16))
            if strcmpi(era_model,'era')
                aps_era_ECMWF_Python(filelist(l,5:16),weatherregstr) %write python donwload file
            elseif strcmpi(era_model,'era5')
                aps_era5_ECMWF_Python(filelist(l,5:16),weatherregstr) %write python donwload file
            end
            python_str = ['python ',filelist(l,5:16),'.py > ',filelist(l,5:16),'down.log &'];
            [a,b] = system(python_str); % start python script
            clear a b
            pause(5);
            cd(workdir)
        else
            if overwrite_flag==-1
                str='';
                while ~strcmpi(str,'y') && ~strcmpi(str,'n') 
                    fprintf(['Do you want to overwrite existing files? \n'])  
                    str = input('[y: for yes, n: no] \n','s');
                end
                if strcmpi(str,'n')
                    overwrite_flag=0;
                else
                    overwrite_flag=1;
                end
            end
            % check if the files need to be overwritten
            if overwrite_flag==1
                cd(subdirpath)
                if exist([subdirpath,'ggap', filelist(l,5:16) '.nc'],'file')~=0
                    delete([subdirpath,'ggap', filelist(l,5:16) '.nc'])
                end
                fprintf('Order and downloading %s \n',filelist(l,5:16))
                if strcmpi(era_model,'era')
                    aps_era_ECMWF_Python(filelist(l,5:16),weatherregstr)    % write python download file
                elseif strcmpi(era_model,'era5')
                    aps_era5_ECMWF_Python(filelist(l,5:16),weatherregstr)   % write python download file
                end
                python_str = ['python ',filelist(l,5:16),'.py > ',filelist(l,5:16),'down.log &'];
                [a,b] = system(python_str); % start python script
                clear a b
                pause(5);
                cd(workdir)
            elseif overwrite_flag==0
                fprintf('File %s has already been downloaded \n',filelist(l,:))
            end
        end
    end
    cd(workdir)
end
cd(workdir)     % back to work directory
    
