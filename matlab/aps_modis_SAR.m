function [] = aps_modis_SAR(batchfile)
% aps_modis_SAR(datelist)
% Scipt to load modis data, mask out clouds, interpolate over gaps and cut 
% the data to the right size for which the tropospheric delay is being computed. 
% The modis data is assumed to be structured in date folders. 
% The batchfile contains the full path to the modis files in these folders. 
% Note that the first line of the batchfile should read "files".
%
% modis calibration factor from:
% Li, Z., Fielding, E. J., Cross, P., & Preusker, R. (2009). 
% Advanced InSAR atmospheric correction: MERIS/MODIS combination and 
% stacked water vapour models. International Journal of Remote Sensing, 30(13), 
% 3343-3363. doi: 10.1080/01431160802562172 
%
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
% OPTIONAL INPUTS:
% conversion            PI conversion factor, default 6.2, or sounding data .
%
% OUTPUTS:
% It will give wet delay maps in cm. 
%
% By David Bekaert - University of Leeds
%
% modifications:
% DB    07/2014     Include calibration factor for modis, dem check
% DB    07/2014     Make more explicit that its only wet delays
% DB    07/2014     Redefine meris_lon(lat)_range to region_lon(lat)_range
% DB    08/2014     Include option to have factors varying for each SAR date
% DB    03/2015     Fix in case file names have varying length
% DB    02/2016     Update an error message




% setting the defaults and checking the input arguments
if nargin<1
    fprintf('aps_modis_SAR(batchfile) \n')
    error('myApp:argChk', ['Not enough input arguments...\n'])
end


conversion_vector = getparm_aps('spectrometer_PIconversion',1);
smpres = getparm_aps('region_res',1);          % in degrees
xlims = getparm_aps('region_lon_range',1);
ylims = getparm_aps('region_lat_range',1);
modis_recalibrated = getparm_aps('modis_recalibrated',1);
stamps_processed = getparm_aps('stamps_processed',1);

% check if the modis data is already calibrated
if strcmpi(modis_recalibrated,'y')
    modis_calibration = 1;
else
    modis_calibration = getparm_aps('modis_calibration',1);
end


%% the actual scripting
%bounds for all ifgms in degrees or in meters
if isempty(xlims) || isempty(ylims)
    error('myApp:argChk', ['Please specify a region_lon_range and region_lat_range\n'])
end
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);
fprintf(['Lon range: ' num2str(xmin) ' -- '  num2str(xmax) ' degrees\n'])
fprintf(['Lat range: ' num2str(ymin) ' -- '  num2str(ymax) ' degrees\n'])
fprintf(['Output resolution is assumed to be ' num2str(smpres) ' degrees \n'])


% getting the number of files to be processed
files = char(textread(batchfile,'%s','headerlines',1));
ix_keep = [];
for counter=1:size(files,1)
    
    if strcmpi(modis_recalibrated,'y')==1
        if ~isempty(strfind(files(counter,:),'recal'));
           ix_keep = [ix_keep; counter];
           recal_str = '_recal';
        end
    else
        if isempty(strfind(files(counter,:),'recal'));
           ix_keep = [ix_keep; counter];
           recal_str = '';

        end
    end
end
files = files(ix_keep,:);
ndates = size(files,1);

% loading the date information  
if strcmp(stamps_processed,'y')
   ps = load(getparm_aps('ll_matfile'));
   ifgs_dates = ps.day;
   fprintf('Stamps processed structure \n')
else
    ifgday_matfile = getparm_aps('ifgday_matfile',1);
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    ifgs_dates = reshape(ifgs_dates,[],1);
    ifgs_dates = unique(ifgs_dates);
end

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

if length(conversion_vector)>1
   fprintf('Conversion factor varies for each SAR date.\n')
end

%start loop here to calculate atmos correction for each date
fprintf('Starting the computation for each SAR date \n')
for n = 1:ndates
    if length(conversion_vector)>1
        ix_date_postion = find(ifgs_dates==datenum(datelist(n,:),'yyyymmdd'));
        conversion =conversion_vector(ix_date_postion);
    else
       conversion =conversion_vector;
    end
    
    file = strtrim([files(n,:)]);
    outfile = [pathlist{n} filesep datelist(n,:) recal_str '_ZWD_nointerp.xyz'];
    outfile_gauss = [pathlist{n} filesep datelist(n,:) recal_str '_ZWD_gauss.xyz'];
    outfile_surf = [pathlist{n} filesep datelist(n,:) recal_str '_ZWD_surf.xyz'];
    [lon_grid,lat_grid,watervap] = grdread2(file);

    D = grdinfo2(file);
    mod_xmin = D(1);
    mod_xmax = D(2);
    mod_ymin = D(3);
    mod_ymax = D(4);
    mod_cols = length(lon_grid);
    mod_rows = length(lat_grid);
    clear D

    %% Calculate wet delay
    %convert from g/cm^2 slant water vapour to cm phase delay 
    watervap=watervap.*conversion.*modis_calibration;
    %mask out dodgy pixels
    grdwrite2(lon_grid,lat_grid,watervap,'tmp.grd')

    % gaussian filter tmp.grd to get tmp2.grd
    grdfil1_cmd = 'grdfilter tmp.grd -Gtmp2.grd -Ni -D2 -Fg50'; %10000
    [a,b]=system(grdfil1_cmd);

    % export tmp.grd as tmp.xyz
    grd2xyz_cmd = ['grd2xyz -R',num2str(mod_xmin),'/',num2str(mod_xmax),'/',num2str(mod_ymin),'/',num2str(mod_ymax),' tmp.grd -bo > tmp.xyz'];
    [a,b]=system(grd2xyz_cmd);

    % use surface to interpolate using tmp.xyz and create tmp_fil.grd
    grdfil_cmd = ['surface -R',num2str(mod_xmin),'/',num2str(mod_xmax),'/',num2str(mod_ymin),'/',num2str(mod_ymax),' tmp.xyz -I',num2str(mod_cols),'+/',num2str(mod_rows), '+ -bi -Gtmp_fil.grd -T0.5'];
    [a,b]=system(grdfil_cmd);

    % downsample and output nointerp, gaussian and surface files
    %nointerp
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp.grd -Gtmp_smp.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo >' outfile];
    [a,b]=system(grd2xyz_cmd);

    %gauss
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp2.grd -Gtmp_smp2.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp2.grd -bo >' outfile_gauss];
    [a,b]=system(grd2xyz_cmd);
    %surf
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp_fil.grd -Gtmp_smp_fil.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp_fil.grd -bo >' outfile_surf];
    [a,b]=system(grd2xyz_cmd);

    % load gaussian and nointerp back in, combine, and output as outfile_gauss
    % opening the data file (not-interpolated)
    nointfid = fopen(outfile,'r');
    data_vector = fread(nointfid,'double');
    fclose(nointfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    noint = data(:,3);
    clear data data_vector


    % opening the gaussian interpolated file
    gaussfid = fopen(outfile_gauss,'r');
    data_vector = fread(gaussfid,'double');
    fclose(gaussfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    gauss = data(:,3);
    xy = data(:,[1:2]);
    clear data data_vector

    % (gaussian interpolated points only where no data available)
    noint(isnan(noint))=0;
    gauss(noint~=0)=0;
    gauss_interp=gauss+noint;
    gauss_interp(gauss_interp==0)=NaN;
    % writing out the date again as a binary table
    data_write = [xy gauss_interp]';
    clear gauss_interp noint gauss
    %output
    fid = fopen(outfile_gauss,'w');
    fwrite(fid,data_write,'double');
    fclose(fid);
    clear data_write

    
    fprintf([num2str(n) ' completed out of ' num2str(ndates) '\n']) 
end

[a,b] = system('rm tmp.grd tmp.xyz tmp2.grd tmp_fil.grd tmp_smp.grd tmp_smp2.grd tmp_smp_fil.grd');

