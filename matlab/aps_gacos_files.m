function [] = aps_gacos_files(orderflag_GACOS_website)
% script that runs and checks which GACOS data files are needed based on
% the satellite pass time.
% INPUTS:
%
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
% By David Bekaert - November 2017
% Modifications:

stdargin = nargin ; 
if stdargin<1
    orderflag_GACOS_website = 0;
end


% getting the variables from the parms_aps file
workdir = pwd;
stamps_processed = getparm_aps('stamps_processed');
UTC_sat =  getparm_aps('UTC_sat');
gacos_datapath = getparm_aps('gacos_datapath');
datestructure = 'yyyymmdd';                               % assumed date structure for era
if isempty(gacos_datapath)
    error('please specify gacos_datapath')
end

% loading the data
if strcmp(stamps_processed,'y')
    ll_matfile = getparm_aps('ll_matfile');
    ps = load(ll_matfile);
    dates = ps.day;
else
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
end

% getting the dates
n_dates = length(dates);

% find two closest times with respect the 1hr spaced ERA5 data
time = str2num(UTC_sat(1:2)) + str2num(UTC_sat(4:5))/60;
fprintf(['Satellite pass is ' num2str(time) ' UTC \n'])
datelist =  datestr(dates,datestructure);
filelist = [datelist repmat('.ztd',size(datelist,1),1) ];

if orderflag_GACOS_website==0
    % outputing this information to a file
    fid = fopen('GACOS_files.txt','w');
    fprintf(['Required GACOS files for all interferograms \n'])
    for k=1:size(filelist,1)
        fprintf([filelist(k,:) '\n']);
        fprintf(fid,[filelist(k,:) '\n']);
    end
    fclose(fid);


    % getting the BBOX needed for GACO, in case user wants to download manually
    crop_range_in = 0.5; % increasing extent of weather data region by this value (degree)
    region_lat_range = getparm_aps('region_lat_range');
    region_lon_range = getparm_aps('region_lon_range');
    if isempty(region_lat_range) == 1
        error('Specify the region for the weather model data')
    end
    fprintf('increasing crop area by %s deg in each direction \n',num2str(crop_range_in))
    S = num2str(min(round(region_lat_range)) - crop_range_in);
    N = num2str(max(round(region_lat_range)) + crop_range_in);
    W = num2str(min(round(region_lon_range)) - crop_range_in);
    E = num2str(max(round(region_lon_range)) + crop_range_in);       % DB fixed typo min to max  
    weatherregstr = [N,'/',W,'/',S,'/',E];   % N/W/S/E
    fprintf('weather model region N/W/S/E %s \n',weatherregstr);


    fid2 = fopen('GACOS_download_info.txt','w');
    % give the metadata needed
    fprintf(fid2,['UTC: ' UTC_sat '\n']);
    fprintf(fid2,['South: ' S '\n']);
    fprintf(fid2,['North: ' N '\n']);
    fprintf(fid2,['West: ' W '\n']);
    fprintf(fid2,['East: ' E '\n\n']);
    fprintf(fid2,['dates:\n']);
    for k=1:size(datelist,1)
        fprintf(fid2,[datelist(k,1:end) '\n']);
    end   
    fclose(fid2);
end


%% Below is specific for the ECMWF data website
crop_range_in = 2; % increasing extent of weather data region by this value (degree)
overwrite_flag=-1;
if orderflag_GACOS_website==1
    
    error('GACOS website does not support api call')

    
    % weather model region
    fprintf('getting the data from the GACOS website ...')
    region_lat_range = getparm_aps('region_lat_range');
    region_lon_range = getparm_aps('region_lon_range');
    if isempty(region_lat_range) == 1
        error('Specify the region for the weather model data')
    end
    fprintf('increasing crop area by %s deg in each direction \n',num2str(crop_range_in))
    S = num2str(min(round(region_lat_range)) - crop_range_in);
    N = num2str(max(round(region_lat_range)) + crop_range_in);
    W = num2str(min(round(region_lon_range)) - crop_range_in);
    E = num2str(max(round(region_lon_range)) + crop_range_in);       % DB fixed typo min to max  
    weatherregstr = [N,'/',W,'/',S,'/',E];   % N/W/S/E
    fprintf('weather model region N/W/S/E %s \n',weatherregstr);

    %%%%% CODE BELOW IS LEGACY CODE FROM ECMWF DOWNLOAD PROGRAM
    %%%%% TB CLEANED
    
    
    
    for l = 1:size(filelist,1)
        subdirpath = [gacos_datapath,'/',filelist(l,5:12),'/'];
        if exist(subdirpath,'dir') == 0
            fprintf('creating directory %s \n',subdirpath);
            mkdir(subdirpath)
        end
        if exist([subdirpath,'ggap', filelist(l,5:16) '.nc'],'file') == 0
            cd(subdirpath)
            fprintf('Order and downloading %s \n',filelist(l,5:16))
            aps_era5_ECMWF_Python(filelist(l,5:16),weatherregstr) %write python donwload file
            python_str = ['python ',filelist(l,5:16),'.py > ',filelist(l,5:16),'down.log &'];
            [a,b] = system(python_str); % start python script
            clear a b
            cd ..
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
                delete([subdirpath,'ggap', filelist(l,5:16) '.nc'])
                fprintf('Order and downloading %s \n',filelist(l,5:16))
                aps_era5_ECMWF_Python(filelist(l,5:16),weatherregstr) %write python donwload file
                python_str = ['python ',filelist(l,5:16),'.py > ',filelist(l,5:16),'down.log &'];
                [a,b] = system(python_str); % start python script
                clear a b
                cd ..
            elseif overwrite_flag==0
                fprintf('File %s has already been downloaded \n',filelist(l,:))
            end
        end
    end
    cd(workdir)
end


% check if the files needs to be organized into a YYYYMMDD folder structure
try
    cd(gacos_datapath)
    for k=1:size(filelist,1)
       file = filelist(k,:);
        date = file(1:8);
        if exist(date,'dir')~=7
            mkdir(date)
        end
        if exist(file,'file')==2
            movefile([file '*'],[date filesep '.'])
        end
        
    end
    cd(workdir)

catch
    cd(workdir)
end


cd(workdir)     % back to work directory
    
