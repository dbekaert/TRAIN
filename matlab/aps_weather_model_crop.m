function [S,N,W,E] = aps_weather_model_crop(crop_range_in)

if nargin<1
    crop_range_in = 2; % increasing extent of weather data region by this value (degree)
end

% weather model region
region_lat_range = getparm_aps('region_lat_range');
region_lon_range = getparm_aps('region_lon_range');
if isempty(region_lat_range) == 1
    error('Specify the region for the weather model data')
end
fprintf('increasing crop area by %s deg in each direction \n',num2str(crop_range_in))
S = num2str(min(round(region_lat_range)) - crop_range_in);
N = num2str(max(round(region_lat_range)) + crop_range_in);
W = num2str(min(round(region_lon_range)) - crop_range_in);
E = num2str(max(round(region_lon_range)) + crop_range_in); 