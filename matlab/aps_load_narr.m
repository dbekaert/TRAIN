function [ Temp,WVapour,Geopot,Pressure,longrid,latgrid,xx,yy,lon0360_flag] = aps_load_narr(file,narr_model,coeff)
% loading narr data, crop to bbox, interpolate onto even grid, and output the variables for TRAIN
% Kyle Murray
%  clear all;close all
load parms_aps
% file = [narr_datapath, filesep, '20150313/20150313_0000'];
load(file);

% Debug figure to test and validate dataloading.
debug_fig = 0;

% Specify the coefficients in case they are not given
if nargin<3
    coeff.Rv= 461.495;           % [j/kg/K] water vapour
    coeff.Rd= 287.05;            % [j/kg/K] specific gas constants dry air
end
if ~isfield(coeff,'Rd') || ~isfield(coeff,'Rv')
    error('Coefficients Rd and Rv are not spcified')
end
%
%% Get the data and interpolate data to even grid

% Define pressure
Plevs =level;

lons=double(lons);
lats=double(lats);
lons2 = lons; 

% Inform about the organisation of the longitudes and add force positive
if sum(lons>180)>1
    lon0360_flag = 'y';
    lonmin=region_lon_range(1);
    lonmax=region_lon_range(2);
else
    lon0360_flag = 'y';
    lons2(lons<0) = lons2(lons<0)+360;
    lonmin=region_lon_range(1)+360;
    lonmax=region_lon_range(2)+360;
end

% Weather model region
crop_range_in = 2; % increasing extent of weather data region by this value (degree)
region_lat_range = getparm_aps('region_lat_range');
region_lon_range = getparm_aps('region_lon_range');
if isempty(region_lat_range) == 1
    error('Specify the region for the weather model data')
end
fprintf('increasing crop area by %s deg in each direction \n',num2str(crop_range_in))
S  = num2str(round(region_lat_range(1)) - crop_range_in);
N = num2str(round(region_lat_range(2)) + crop_range_in);
W = num2str(round(lonmin) - crop_range_in);
E = num2str(round(lonmax) + crop_range_in);       % DB fixed typo min to max
weatherregstr = [N,'/',W,'/',S,'/',E];   % N/W/S/E
fprintf('weather model region: N/W/S/E %s \n',weatherregstr);

% Make all longitudes have positive values
x = str2num(W):0.25:str2num(E);
y = str2num(S):0.25:str2num(N);    % Create a 1/4 degree lat-lon grid

[xx1,yy1]=meshgrid(x,y);
n_p=length(air(:,1,1))*length(air(1,:,1));
% air_intrp = griddata(double(lons2(:)),double(lats(:)),reshape(air(:,:,ii),[n_p,1]),xx,yy);  % interpol
% figure;imagesc(x,y,air_intrp);colorbar;axis image

for ii=1:length(Plevs) %loop through pressures
    Temp(:,:,ii) = griddata(double(lons2(:)),double(lats(:)),reshape(air(:,:,ii),[n_p,1]),xx1,yy1);  % interpol
    qv(:,:,ii) = griddata(double(lons2(:)),double(lats(:)),reshape(shum(:,:,ii),[n_p,1]),xx1,yy1);  % interpol
    H(:,:,ii) = griddata(double(lons2(:)),double(lats(:)),reshape(hgt(:,:,ii),[n_p,1]),xx1,yy1);  % interpol
end
% Convert Geopotential Height to geopotential
%  important this is needed to be consitent with the other weather modules
g0 = 9.80665;
Geopot = H.*g0;

% Update the no data values with NaN's
Temp(Temp==missing_value)=NaN;
qv(qv==missing_value)=NaN;
Geopot(Geopot==missing_value)=NaN;

% Number of lon lat grid nodes
n_latitude_points = length(yy1(:,1));
n_longitude_points = length(yy1(1,:));

% Generate the pressure levels
Pressure = repmat(Plevs,[1,n_longitude_points,n_latitude_points]);
Pressure = permute(Pressure,[3,2,1]);

% Computation of the saturated water vapour e
% use mixing ratio r = E*e/(p-e), where E = Rd/Rv
% and specific humidity qv = r/(1+r)
E = coeff.Rd./coeff.Rv;
WVapour= qv.*Pressure./(E.*(1-qv)+qv);
WVapour = reshape(WVapour,size(qv));
clear qv

% Get list of points to look at in analysis
[xx,yy] = meshgrid(1:n_longitude_points,1:n_latitude_points);

% generate the lon lat grid
latgrid = repmat(yy1,[1,1,length(Plevs)]);
longrid = repmat(xx1,[1,1,length(Plevs)]);

%% validation plots
if debug_fig ==1
    figure('position',[ 3         628        1402         586]);
    subplot(2,5,1)
    imagesc(Temp(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('temp upper atmo')
    subplot(2,5,2)
    imagesc(Temp(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('temp lower atmo')
    
    subplot(2,5,3)
    imagesc(Pressure(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('pressure upper atmo')
    subplot(2,5,4)
    imagesc(Pressure(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('pressure lower atmo')
    
    
    
    subplot(2,5,6)
    imagesc(WVapour(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('water vapour upper atmo')
    subplot(2,5,7)
    imagesc(WVapour(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('water vapour lower atmo')
    
    
    subplot(2,5,8)
    imagesc(Geopot(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('geopotential upper atmo')
    subplot(2,5,9)
    imagesc(Geopot(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('geopotential lower atmo')
    
    
    subplot(2,5,5)
    imagesc(latgrid(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('lat')
    subplot(2,5,10)
    imagesc(longrid(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('lon')
end
