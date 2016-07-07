function [] = aps_spectrometer_PWV_modis(batchfile)
% aps_spectrometer_PWV_modis(batchfile)
% Scipt to load modis data, mask out clouds. The modis data is assumed to 
% be structured in date folders. The batchfile contains the full path to the
% modis files in these folders. Note that the first line of the batchfile 
% should read "files".
%

% INPUTS used:
% batchfile             A txt file containing the full path and file names of the
%                       meris data that needs to be processed. The first
%                       line of this file should read "files". The data
%                       should be structured in date folders.
% xlims                 Limits in the x-direction, either in degrees
% ylims                 Limits in the y-direction, either in degrees
% smpres                The output resolution, either in degrees
%                       Units needs to be consistend with xlims and ylims.
%
% By David Bekaert - University of Leeds
% August 2014

% setting the defaults and checking the input arguments
if nargin<1
    fprintf('aps_spectrometer_PWV_modis(batchfile) \n')
    error('myApp:argChk', ['Not enough input arguments...\n'])
end


smpres = getparm_aps('region_res',1);          % in degrees
xlims = getparm_aps('region_lon_range',1);
ylims = getparm_aps('region_lat_range',1);
stamps_processed = getparm_aps('stamps_processed',1);


%% the actual scripting
%bounds for all ifgms in degrees or in meters
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);


% getting the number of files to be processed
files = char(textread(batchfile,'%s','headerlines',1));
ix_keep = [];
for counter=1:size(files,1)
    if isempty(strfind(files(counter,:),'recal'));
       ix_keep = [ix_keep; counter];
    end
end
files = files(ix_keep,:);
% exclude recalibrated files
ndates = size(files,1);

% loading the date information  
if strcmp(stamps_processed,'y')
   ps = load(getparm_aps('ll_matfile',1));
   ifgs_dates = ps.day;
   fprintf('Stamps processed structure \n')
else
    ifgday_matfile = getparm_aps('ifgday_matfile',1);
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    ifgs_dates = reshape(ifgs_dates,[],1);
    ifgs_dates = unique(ifgs_dates);
end

fprintf(['Lon range: ' num2str(xmin) ' -- '  num2str(xmax) ' degrees\n'])
fprintf(['Lat range: ' num2str(ymin) ' -- '  num2str(ymax) ' degrees\n'])
fprintf(['Output resolution is assumed to be ' num2str(smpres) ' degrees \n'])

% extracting the dates from the filenames
for k=1:ndates
    [path,filename_temp,ext_temp] = fileparts(files(k,:));
    clear filename_temp ext_temp
    [path_temp,date,ext_temp] = fileparts(path);
    clear path_temp ext_temp
    
    % save the paths as structures to allow for variable path lengths
    pathlist{k} = path;
    clear path
    % saving the date information
    datelist(k,:) =date;
    clear date
end

%start loop here to calculate atmos correction for each date
fprintf('Starting the masking and writing of the PWV for each SAR date \n')
for n = 1:ndates
     
    file = [files(n,:)];
    outfile_watervapor = [pathlist{n} filesep datelist(n,:) '_ZPWV_nointerp.xyz'];
    [lon_grid,lat_grid,watervap] = grdread2(file);          % g/cm^2 water vapour

    D = grdinfo2(file);
    mod_xmin = D(1);
    mod_xmax = D(2);
    mod_ymin = D(3);
    mod_ymax = D(4);
    mod_cols = length(lon_grid);
    mod_rows = length(lat_grid);
    clear D

    %% Calculate wet delay
    %mask out dodgy pixels
    grdwrite2(lon_grid,lat_grid,watervap,'tmp.grd')

    % downsample and output without interpolation
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' -F tmp.grd -Gtmp_smp.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo >' outfile_watervapor];
    [a,b]=system(grd2xyz_cmd);


    % load nointerp back in and write out
    nointfid = fopen(outfile_watervapor,'r');
    data_vector = fread(nointfid,'double');
    fclose(nointfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    noint = data(:,3);
    xy = data(:,[1:2]);
    clear data data_vector

%     figure; scatter3(xy(:,1),xy(:,2),noint,15,noint,'filled'); view(0,90); axis equal; axis tight

    %output
    data_write = [xy noint]';
    fid = fopen(outfile_watervapor,'w');
    fwrite(fid,data_write,'double');
    fclose(fid);
    clear data_write


    
    fprintf([num2str(n) ' completed out of ' num2str(ndates) '\n']) 
end

[a,b] = system('rm tmp.grd tmp.xyz tmp2.grd tmp_fil.grd tmp_smp.grd tmp_smp2.grd tmp_smp_fil.grd');

