function [] = aps_wrf(start_step,end_step)
% [] = aps_wrf(start_step,end_step)
% Function that computes the inteferometric tropopsheric delays from
% WRF simulations
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
% 11/2013 	DB 	Include script to setup WRF
% 11/2013 	DB 	Include script to make WPS namelist.wps files
% 04/2015   DB  Remove test version of aps_wrf_SAR.m with stable version


% current processing directory
curdir = pwd;


if start_step==0
    % Setting up the WRF input files needed
    fprintf('Step 0: Setting up the WPS namelist.wps and WRF namelist.input files for each SAR date, \n compiles list of required datafiles (WRF_ files.txt). \n')
    aps_wrf_files
end
if start_step==0.1
    fprintf('Setting up the WRF filestructure \n')
    aps_wrf_files_setup
end
if start_step==1
    % SAR delays using WRF data
    fprintf('Step 1: Compute the WRF tropospheric (zenith) delay for individual dates\n')
    % running the WRF computation
    aps_wrf_SAR

end

if start_step<=2 && end_step >=2
    fprintf('Step 2: Computes WRF tropospheric (slant) delay for inteferograms\n')
    aps_wrf_InSAR   
end

cd(curdir)


