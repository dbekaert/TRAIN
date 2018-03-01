function [ Temp,e,Geopot,P,longrid,latgrid] = aps_weather_model_nan_check( Temp,e,Geopot,P,longrid,latgrid);
% Function which fills in the gaps in the weather model data with the
% nearest neighbor value, or from the level above.
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
% April 2016


% getting the pressure 
Pressure_level = squeeze(P(1,1,:));
lon = double(longrid(:,:,1));
lat = double(latgrid(:,:,1));

% step_levels
% define it in such a way we loop from the upper atmopshere to the lower atmosphere
step_level = [length(Pressure_level):-1:1]';
if Pressure_level(1)<Pressure_level(end) % first level is top APS
    step_level = [1:length(Pressure_level)]';
end

% looping though and use TEMP as inital check for nan
n_pixels = size(P,1)*size(P,2);
for k=1:length(Pressure_level)
    ix_nan = isnan(Temp(:,:,step_level(k)));
    n_nan = sum(sum(ix_nan));
    
    % check if there are any NaN values
    if n_nan<n_pixels-3
        temp_data = Temp(:,:,step_level(k));
        Temp(:,:,step_level(k)) = griddata(reshape(lon(~ix_nan),[],1),reshape(lat(~ix_nan),[],1),double(reshape(temp_data(~ix_nan),[],1)),lon,lat,'nearest'); 
   
        temp_data = e(:,:,step_level(k));
        e(:,:,step_level(k)) = griddata(reshape(lon(~ix_nan),[],1),reshape(lat(~ix_nan),[],1),double(reshape(temp_data(~ix_nan),[],1)),lon,lat,'nearest'); 
        
    % check if all are NaN values
    elseif n_nan==n_pixels
        if k==1
            error('Weird seems like the top of the atmopshere has no data')
        end
        Temp(:,:,step_level(k))=Temp(:,:,step_level(k-1));
        e(:,:,step_level(k))=e(:,:,step_level(k-1));
    else
        Temp(:,:,step_level(k))=nanmean(nanmean(Temp(:,:,step_level(k))));
        e(:,:,step_level(k))=nanmean(nanmean(e(:,:,step_level(k))));
        
    end
end


