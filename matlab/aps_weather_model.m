function [] = aps_weather_model(model_type,start_step,end_step,save_path)
% [] = aps_weather_model(model_type,start_step,end_step)
% Function that computes the interferometric tropospheric delays from
% weather model including ERA-Interim, ERA5 (test data), NARR, MERRA, MERRA2.
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
% November 2013
%
% Modifications:
% 02/2014   DB  Include data from ECMWF and BADC as option
% 04/2014   DB  Include auto download for ECMWF website
% 04/2016   DB  Convert to a generic weather model correction using modular
%               approach. Include support for MERRA.
% 07/2017   DB 	Adding ERA5 model based on test-data
% 02/2018   KM 	Adding NARR model based on test-data

% current processing directory
curdir = pwd;

% error catching
if nargin <3
    error('Need to specify at least model_type, start_step, end_step');
end
if ~strcmpi(model_type,'era5') & ~strcmpi(model_type,'era') & ~strcmpi(model_type,'gacos') & ~strcmpi(model_type,'merra') & ~strcmpi(model_type,'narr') & ~strcmpi(model_type,'merra2')
    error(['model_type needs to be era, merra, merra2, gacos, era5, narr'])
end
% define save path if not given
if nargin<4
   save_path = [curdir filesep]; 
end


%% Defining the save path in case not given
% filenames
curdir = pwd;
model_type = lower(model_type);
% saving part of the data in a subfolder
save_path = [save_path filesep 'aps_' model_type];
if exist(save_path,'dir')~=7
    mkdir(save_path);
end

%% The different steps in the code
if start_step==0
    if strcmpi(model_type,'era') 
        % Dummy run on the needed ERA-I files
        fprintf('Step 0: Dummy run on the needed ERA-I data files \n')
        aps_era_files(0);
     elseif strcmpi(model_type,'era5')
        % Dummy run on the needed ERA5 files
        fprintf('Step 0: Dummy run on the needed ERA5 data files \n')
        aps_era5_files(0);
     elseif strcmpi(model_type,'gacos')
        % Dummy run on the needed gacos files
        fprintf('Step 0: Dummy run on the needed gacos data files \n')
        aps_gacos_files(0);
     elseif strcmpi(model_type,'merra') || strcmpi(model_type,'merra2')
        % required MERRA files
        aps_merra_files(0,model_type);
     elseif strcmpi(model_type,'narr') 
        % required narr files
        aps_narr_files(0);
    end
end
if start_step<=1 && end_step >=1 
    fprintf(['Step 1: Order and Download ' upper(model_type) ' data files \n'])
    if strcmpi(model_type,'era') 
        % order the ECMWF data ERA-I files
        era_data_type = getparm_aps('era_data_type',1);
        if strcmpi(era_data_type,'ECMWF')
            aps_era_files(1);
        else
            fprintf('You need to download BADC data from command line')
        end
    elseif strcmpi(model_type,'gacos')
        fprintf('You need to download GACOS data from the website, use GACOS_download_info.txt for relevant information\n')
    elseif strcmpi(model_type,'era5')
        aps_era5_files(1);
    elseif strcmpi(model_type,'merra') || strcmpi(model_type,'merra2')
        aps_merra_files(1,model_type);
    elseif strcmpi(model_type,'narr') 
        aps_narr_files(1);
    end
end
if start_step<=2 && end_step >=2 
    % SAR zenith delays 
    fprintf(['Step 2: Compute the ' upper(model_type) ' tropospheric (zenith) delay for individual dates\n'])
    % running the SAR delay computation using the weather model observations 
    if strcmpi(model_type,'gacos')
        fprintf('GACOS delays as downloaded are already Zenith delays, do not need to do anything here\n')
    else
        aps_weather_model_SAR(model_type);
    end
end

if start_step<=3 && end_step >=3
    % if GACOS check if the files are in date structure
    if strcmpi(model_type,'gacos')
        aps_gacos_files(-1);
    end
    
    % InSAR slant delays
    fprintf(['Step 3: Computes '  upper(model_type) ' tropospheric (slant) delay for inteferograms\n'])
    aps_weather_model_InSAR(model_type);
end

cd(curdir)


