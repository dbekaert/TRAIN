function [modelfile_before,modelfile_after] = aps_weather_model_filenames(model_type,time_before,time_after,date_before, date_after,weather_model_datapath)
% function which generates the filenames of the weather model data
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
% inputs:
% model_type either era or merra
% time_before, time_after are strings 'HHMM' with the weather model times before and after UTC acquisition time.
% date_before, date_after are the dates 'YYYYMMDD' with the weather model data before and after SAR acqusition
%
% modification:
% DB    02/05/2016  Fix bug for MERRA and MERRA 2, as the before and after were defined the same
% DB    12/11/2017  Update merra to be netcdf 4, hdf5 does not seem to work
%                   after ftp server changed.
% KM    02/26/2018  Including NARR support

if nargin<6
    weather_model_datapath = pwd;
end

for d =1:size(date_before,1)
    if strcmpi(model_type,'era')    
        %Format ggapYYYYMMDDHHMM.nc
        modelfile_before(d,:) = [weather_model_datapath filesep date_before(d,:) filesep 'ggap' date_before(d,:) time_before(d,:) '.nc']; 
        modelfile_after(d,:) = [weather_model_datapath filesep date_after(d,:) filesep 'ggap' date_after(d,:) time_after(d,:) '.nc'];
    elseif strcmpi(model_type,'merra') 
        %Format MERRA_YYYYMMDD_HH.hdf
        %modelfile_before(d,:) = [weather_model_datapath filesep date_before(d,:) filesep 'MERRA_' date_before(d,:) '_' time_before(d,1:2) '.hdf'];         
        %modelfile_after(d,:) = [weather_model_datapath filesep date_after(d,:) filesep 'MERRA_' date_after(d,:) '_' time_after(d,1:2) '.hdf'];
        modelfile_before(d,:) = [weather_model_datapath filesep date_before(d,:) filesep 'MERRA_' date_before(d,:) '_' time_before(d,1:2) '.nc4'];         
        modelfile_after(d,:) = [weather_model_datapath filesep date_after(d,:) filesep 'MERRA_' date_after(d,:) '_' time_after(d,1:2) '.nc4'];
    elseif strcmpi(model_type,'merra2')
        %Format MERRA2_YYYYMMDD_HH.nc4
        modelfile_before(d,:) = [weather_model_datapath filesep date_before(d,:) filesep 'MERRA2_' date_before(d,:) '_' time_before(d,1:2) '.nc4'];         
        modelfile_after(d,:) = [weather_model_datapath filesep date_after(d,:) filesep 'MERRA2_' date_after(d,:) '_' time_after(d,1:2) '.nc4'];
    elseif strcmpi(model_type,'narr')
        %Format YYYYMMDD_HHMM.mat
        modelfile_before(d,:) = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_' time_before(d,1:4) '.mat'];         
        modelfile_after(d,:) = [weather_model_datapath filesep date_after(d,:) filesep date_after(d,:) '_' time_after(d,1:4) '.mat'];

    end
end