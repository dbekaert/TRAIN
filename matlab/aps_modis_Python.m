function []=aps_modis_Python()
% function that downloads the MOSID data from OSCAR and saves it in the
% correct data structure
%
% A loop is included rather than a simultaneous call as this would kill the
% OSCAR server. To make sure the service was succesfull, a check is
% performed till all files are downloaded.
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
% By David Bekaert - University of Leeds - May 2014 
%
% modifications
% DB    07/2014     Redefine meris_lat(lon)_range to region_lat(lon)_range



% getting the variables from the parms_aps file
workdir = pwd;
stamps_processed = getparm_aps('stamps_processed');
UTC_sat =  getparm_aps('UTC_sat');
modis_datapath = getparm_aps('modis_datapath');
datestructure = 'yyyymmdd';                               % assumed date structure for era

% loading the data
if strcmp(stamps_processed,'y')
    ll_matfile = getparm_aps('ll_matfile');
    ps = load(ll_matfile);
    dates = ps.day;
    load psver
else
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
end    

% the region the delay needs to cover.
% slightly larger than the InSAR region
xlims = getparm_aps('region_lon_range');
ylims = getparm_aps('region_lat_range');
smpres = getparm_aps('region_res');
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);


% the region which is cropped from the ERA data and used to make the interpolation.
% Should be  larger than the region to which the delay is computed
lonmin = floor(xmin)-1;
lonmax= ceil(xmax)+1;
latmin = floor(ymin)-1;
latmax = ceil(ymax)+1;



% getting the dates
n_dates = length(dates);

% downloading of the data
counter=1;
continueflag=1;
count_exist=0;
while continueflag==1


    for k=1:n_dates
        if exist([modis_datapath filesep datestr(dates(k),datestructure)],'dir')~=7
            mkdir([modis_datapath filesep datestr(dates(k),datestructure)]);
        end
        if exist([modis_datapath filesep datestr(dates(k),datestructure) filesep 'OSCAR_Modis_' datestr(dates(k),datestructure) '.grd'  ],'file') ~=2

            cd([modis_datapath filesep datestr(dates(k),datestructure)]);

            fprintf(['Downloading modis data for ' datestr(dates(k),datestructure) '\n'])
            
            %  options for get_modis
            % -r (--region): lon in [-180,180], lat in [-90,90], mandatory'
            % -t (--time): time and date in ISO format, mandatory'
            % -p (--platform): terra (default), aqua, or any.'
            % -o (--outfile): output filename (default: mod_minlon_maxlon_minlat_maxlat_dateTtime.grd)'
            % -w (--timewindow): saerch time window size in second (default: 18000, i.e. time +/-5hrs)'
            % -g (--gridsize): grid spacing in degrees (default: 30./3600.=0.008333333...)'
            % -l (--localtime): if used, the given time is local time, otherwise UTC'
            % -f (--figurefile): downloads the png figure file (same fileroot but with .png)'
            % -v (--verbose)'
            % -s (--server): specify which server to use - oscar1 (default) or oscar2'

            commandstr = ['python $get_modis_filepath -r ' num2str(lonmin) '/' num2str(lonmax) '/' num2str(latmin) '/' num2str(latmax) ' -t ' datestr(dates(k),'yyyy-mm-dd') 'T' UTC_sat ':00 -p terra -v -o OSCAR_Modis_' datestr(dates(k),datestructure) '.grd -f -v > download.log' ]; 
            [a,b] = system(commandstr);
            clear commandstr a b
         else
            count_exist = count_exist+1;
         end
    end

    if count_exist==n_dates
        fprintf('All files have been downloaded \n')
        continueflag=0;
    end
    
    if counter==10
        fprintf('Stop iterating the download. Not all SAR dates where downloaded... \n')
        continueflag=0;
    end

counter = counter+1;
end
cd(workdir)




   

