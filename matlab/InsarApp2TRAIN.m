% Set-up script for ISCE processing to TRAIN.
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
% By David Bekaert 
% modifications:
% DB    20/12/2012      Include loading and storing of the IFG's too for
%                       phase-based methods

close all
clear all
clc

%% SETUP SCRIPT compatible with InsarApp.py
% inputs from user
% IFG folders are "date1-date2" or "date1_date2". 
% How are date1 and date2 defined?
date_format = 'yyyymmdd';           % this needs to be recognized by matlab.
                                    % mmm for jan format, MM is min two digit, 
                                    % while hh is hour 2 digit. You can
                                    % also include - or / if needed. e.g.
                                    % 'yyyy-mmm-dd'
                                    
% MASK file. In case there is a file called "MASK" with the same
% size of the data then it will be applied. 
% in addtion, nodata in the lonlat data can be masked as well.
nodata_lonlat =0;                   % either a number of [] to ignore

% also save the phase data in a stamps-like format which allows to use plotting
% functionality from stamps.
save_data = 'y';                    % string either 'y' or 'n' 

% radar or geo-coordinates
process_type = 'geo';               % a string either 'radar' or 'geo'



% satellite UTC time
UTC_time = '14:30';                 % a string defined as 'HH:MM'

% full path to your DEM file. Easiest is if this is a grd. 
DEMfile = '/Volumes/Bekaert_disk2/Projects2/DEM/Taiwan/dem.grd';

% TRAIN directory (best not to change such it remains consistent between programs)
TRAIN_dir = 'aps_train';            % the location where the TRAIN dir will be set-up

% plot debug figures
plot_debug = 'y';                   % string either 'y' or 'n' 

% filename of the unwrapped data
unw_filename = 'filt_topophase.unw';    % note the '.geo' string is added automatically if needed

%% CHECK if the DEM file exist
if exist(DEMfile,'file')~=2
    error([DEMfile ' does not exist']);
end

%% check if the type is correct
geo_str = '';
if strcmpi(process_type,'geo')
    geo_str = '.geo';
end
if ~strcmpi(process_type,'geo') && ~strcmpi(process_type,'radar') 
    error(['process_type can only be ''geo'' or ''radar'''])
end


%% check on the filenames if needed
[temp1,temp2,temp3] =fileparts(unw_filename);
if strcmpi(geo_str,'.geo') && ~strcmpi(temp3,'.geo')
    fprintf(['Looks like you want to run train in .geo mode\nYour unw_filename does not have a ".geo" extension\n'])
    fprintf(['***You can ignore this if you do not plan to run the TRAIN phase-based correction methods.***\n'])
    % Check if action is needed
    repeat=1;
    while repeat==1
        action_flag= input('Do you want to update it to .geo? [y/n] ','s');
        if strcmpi(action_flag,'y')
            repeat=0;
        elseif strcmpi(action_flag,'n')
            repeat=0;
        end
    end
    if strcmpi(action_flag,'y')
        unw_filename = [unw_filename '.geo'];
    else
        repeat=1;
        while repeat==1
            action_flag= input('Do you want to continue? [y/n] ','s');
            if strcmpi(action_flag,'y')
                repeat=0;
            elseif strcmpi(action_flag,'n')
                repeat=0;
            end
        end
        if strcmpi(action_flag,'n')
           error('Abord by user. Update the unw_filename') 
        end
    end
end

%% generate output directory to save along the way data
TRAIN_dir = [pwd filesep TRAIN_dir filesep 'SMALL_BASELINES' filesep];
if exist(TRAIN_dir,'dir')~=7
    mkdir(TRAIN_dir)
else
    error(['Directory exist: remove and rerun: ' TRAIN_dir])
end

%% GENERATE a list of IFG COMBINATIONS
list1 = dir('1*');      % get all the dir pre 2000
list2 = dir('2*');      % get all the dir post 2000
list = [list1 ;list2];  % combine the lists
clear list1 list2

% location of the lon and lat file
latfile = [];
latxmlfile = [];
lonfile = [];
lonxmlfile = [];
losfile = [];
losxmlfile = [];
hgtfile = [];
hgtxmlfile = [];
unwfile = [];
% information on the sensor
insarappxmlfile = [];
% extracting the ifg combinations
for k = 1:length(list)
    ifgday(k,:) = list(k).name;
    
    % search for a lon, lat and xml file
    if exist([ ifgday(k,:) filesep 'lon.rdr' geo_str],'file')==2 && exist([ ifgday(k,:) filesep 'lat.rdr' geo_str],'file')==2 && exist([ ifgday(k,:) filesep 'lon.rdr' geo_str '.xml'],'file')==2 && exist([ ifgday(k,:) filesep 'lat.rdr' geo_str '.xml'],'file')==2 && isempty(lonfile)
        lonfile = [ ifgday(k,:) filesep 'lon.rdr' geo_str];  
        lonxmlfile = [ ifgday(k,:) filesep 'lon.rdr' geo_str '.xml'];
        latfile = [ ifgday(k,:) filesep 'lat.rdr' geo_str];
        latxmlfile = [ ifgday(k,:) filesep 'lat.rdr' geo_str '.xml'];
    end
    % search for the los file
    if exist([ ifgday(k,:) filesep 'los.rdr' geo_str],'file')==2 && exist([ ifgday(k,:) filesep 'los.rdr' geo_str '.xml'],'file')==2 && isempty(losfile)
        losfile = [ ifgday(k,:) filesep 'los.rdr' geo_str];  
        losxmlfile = [ ifgday(k,:) filesep 'los.rdr' geo_str '.xml'];
    end
    % search for the hgt file
    if exist([ ifgday(k,:) filesep 'z.rdr' geo_str],'file')==2 && exist([ ifgday(k,:) filesep 'z.rdr' geo_str '.xml'],'file')==2 && isempty(hgtfile)
        hgtfile = [ ifgday(k,:) filesep 'z.rdr' geo_str];  
        hgtxmlfile = [ ifgday(k,:) filesep 'z.rdr' geo_str '.xml'];
    end
    % search for the InSAR app xml file
    if exist([ ifgday(k,:) filesep 'insarApp.xml'],'file')==2 && isempty(insarappxmlfile)
        insarappxmlfile = [ ifgday(k,:) filesep 'insarApp.xml'];        
    end    
    % search for the unwrapped interferogram files
    if exist([ ifgday(k,:) filesep unw_filename],'file')==2
        unwfile{k} = [ifgday(k,:) filesep unw_filename];
    else
        unwfile{k} = [];
    end
end
clear list
% removing the separation between folders
if length(ifgday)<2*length(date_format)
    error('You folder names are not of format "DATE1[-/_/ /]DATE2')
end
ifgday1 = ifgday(:,1:length(date_format));
ifgday2 = ifgday(:,end-length(date_format)+1:end);
clear ifgday

% converting the dates in numerical dates
ifgday = [str2num(datestr(datenum(ifgday1,date_format),'yyyymmdd')) str2num(datestr(datenum(ifgday2,date_format),'yyyymmdd')) ];
clear ifgday1 ifgday2
ifgday = [datenum(num2str(ifgday(:,1)),'yyyymmdd') datenum(num2str(ifgday(:,2)),'yyyymmdd')];
day = unique(ifgday);
n_ifg = size(ifgday,1);
ifgday_ix = [];
for k=1:n_ifg
   ifgday_ix(k,1) = find(ifgday(k,1)==day);
   ifgday_ix(k,2) = find(ifgday(k,2)==day);
end
master_ix = 1;
master_day = day(master_ix);

%% GET LONLAT coordinates 
if isempty(lonfile) || isempty(latfile) || isempty(latxmlfile)
    error('Could not file lon.rdr, lat.rdr, or lat.rdr.xml file')
end

% getting the file WIDTH's
[latwidth] = get_parm_xml(latxmlfile,'width');
[lonwidth] = get_parm_xml(lonxmlfile,'width');
[loswidth] = get_parm_xml(losxmlfile,'width');
if latwidth~=lonwidth || latwidth~=loswidth
    error(['Width of the lon, lat and los are not consistent'])
else
    WIDTH = latwidth;
    clear latwidth lonwidth lonwidth;
end
fprintf(['Width is ' num2str(WIDTH) '\n'])

% loading the lonlat information
% loading lon
LON = load_isce(lonfile);
if strcmpi(plot_debug,'y')
    figure('name','LON');
    imagesc(LON)
    axis equal; axis tight; colorbar; axis xy
end
    
% loading lat
LAT = load_isce(latfile);
if strcmpi(plot_debug,'y')
    figure('name','LAT');
    imagesc(LAT)
    axis equal; axis tight; colorbar; axis xy
end

% check if the lonlat needs to be masked
MASKED_flag = 'n';
if ~isempty(nodata_lonlat)
    MASK = ~((LON==nodata_lonlat) | (LAT==nodata_lonlat));        % mask of points to keep;
    MASKED_flag = 'y';
else
    MASK = logical(ones(size(LON)));
end

% check if a MASK file exist
if exist('MASK','file')==2
    MASK_temp = load_isce(MASK); 
    fprintf('Need to write piece of code which combines these MASK files together \n ')
    MASKED_flag = 'y';

    keyboard
end

% combining
lonlat = [LON(MASK) LAT(MASK)];
%lonlat = [reshape(LON,[],1) reshape(LAT,[],1)];
n_ps = size(lonlat,1);
clear LON LAT


% getting a position file that relates a vector to matrix position
[LENGTH] = get_parm_xml(latxmlfile,'length');
[J,I] = meshgrid([1:LENGTH],[1:WIDTH]);
IND = sub2ind([WIDTH LENGTH], I(MASK), J(MASK));
%IND = sub2ind([WIDTH LENGTH], reshape(I,[],1), reshape(J,[],1));
clear I J 
% generate template xml for aps file
templatexmlfile = lonxmlfile;
aps_save([TRAIN_dir 'isce2train.mat'],IND,geo_str,WIDTH,LENGTH,MASKED_flag,templatexmlfile)


%% getting the heights if found
phase_based_good = 0;
if ~isempty(hgtfile)
    fprintf('Found heights and will store...\n')
    HGT = load_isce(hgtfile); 
    hgt = HGT(MASK);
    clear HGT;
    aps_save([TRAIN_dir 'hgt2.mat'],hgt)
    phase_based_good = 1;
end

%% getting the unwrapped interferograms if found
if ~isempty(unwfile)
    fprintf('Found unwrapped IFGs, will store...\n')
    % try and allocate the memory in advance
    ph_uw = zeros([n_ps n_ifg]);
    for k=1:length(unwfile)
        if exist(unwfile{k},'file')==2
            try
                [DATA] = load_isce(unwfile{k});
                DATA = DATA(:,:,2);                 % first band is amplitude, second band is unwrapped phase
                % apply the mask and store the data
                ph_uw(:,k) = DATA(MASK);
                clear DATA
            catch
                keyboard
            end
        end
    end
    phase_based_good = phase_based_good+1;
    aps_save([TRAIN_dir 'phuw_sb2.mat'],ph_uw)
end

if phase_based_good==2
   fprintf('Found all the information to run phase-based TRAIN corrections\n') 
   phase_based_good = 'y';
else
    phase_based_good = 'n';
end
aps_save([TRAIN_dir 'isce2train.mat'],phase_based_good)

%% getting the incidence angle
[TEMP] = load_isce(losfile);
INC = TEMP(:,:,1);
HEADING = TEMP(:,:,2);
clear TEMP
if strcmpi(plot_debug,'y')
    figure('name','INC');
    imagesc(INC)
    axis equal; axis tight; colorbar; axis xy
    figure('name','HEADING');
    imagesc(HEADING)
    axis equal; axis tight; colorbar; axis xy
end
la = INC(MASK)./180*pi;
%la = reshape(INC,[],1)./180*pi;
clear INC
heading = nanmean(HEADING(MASK));
clear HEADING

%% getting the sensor information 
% getting the file WIDTH
[SENSOR] = get_parm_xml(insarappxmlfile,'Sensor name');
fprintf(['Sensor is ' num2str(SENSOR) '\n'])

% getting the wavelength hard-coded
if strcmpi(SENSOR,'alos')
    lambda = 23.6;
elseif strcmpi(SENSOR,'radarsat2') || strcmpi(SENSOR,'radarsat') || strcmpi(SENSOR,'rsat') || strcmpi(SENSOR,'rsat2') ||   strcmpi(SENSOR,'envisat') || strcmpi(SENSOR,'ers') || strcmpi(SENSOR,'env') || strcmpi(SENSOR,'ers1')  || strcmpi(SENSOR,'ers2')
    lambda = 5.6;
else
   error('Sensor wavelength not yet hard-coded') 
end
% convert to m
lambda =lambda./100;



%% saving the data
% stamps-like format
psver =2;
save([TRAIN_dir 'psver.mat'],'psver')
aps_save([TRAIN_dir 'ps2.mat'],lonlat,n_ifg,n_ps,day,ifgday_ix,ifgday,master_day,master_ix)
aps_save([TRAIN_dir 'la2.mat'],la)
dlmwrite([TRAIN_dir 'lambda.1.in'],lambda)
dlmwrite([TRAIN_dir 'heading.1.in'],heading)
dlmwrite([TRAIN_dir 'WIDTH.txt'],WIDTH)
fid = fopen([TRAIN_dir 'processor.txt'],'w');
fprintf(fid,'isce');
fclose(fid)
% will make few dummy directories
ERA_dir = [TRAIN_dir 'ERA'];
MERRA_dir = [TRAIN_dir 'MERRA'];
MODIS_dir = [TRAIN_dir 'MODIS'];
MERIS_dir = [TRAIN_dir 'MERIS']; 
if exist([ERA_dir],'dir')~=7
    mkdir(ERA_dir)
end
if exist([MERRA_dir],'dir')~=7
    mkdir(MERRA_dir)
end
if exist([MERIS_dir],'dir')~=7
    mkdir(MERIS_dir)
end
if exist([MODIS_dir],'dir')~=7
    mkdir(MODIS_dir)
end
cd(TRAIN_dir)
% setting the stamps parameters
getparm                                 % this is for PS
setparm('small_baseline_flag','y')      % change to SB    
setparm('percent_rand',[],-1)           % resetting all PS to SB variables
setparm('weed_standard_dev',[],-1)
setparm('scla_drop_index',[],-1)
setparm('unwrap_method',[],-1)
setparm('merge_resample_size',[],-1)
getparm                                 % now this will be for SB
% setting the TRAIN parameters
getparm_aps
setparm_aps('UTC',UTC_time)
setparm_aps('era_datap',ERA_dir)
setparm_aps('merra_datap',MERRA_dir)
setparm_aps('meris_datap',MERIS_dir)    
setparm_aps('modis_datap',MODIS_dir)    
setparm_aps('demfile',DEMfile)    




