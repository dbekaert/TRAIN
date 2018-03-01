function [timelist,model_step] = aps_weather_model_time(model)
% script that return the model time spacing
% INPUTS:
% 
%     Copyright (C) 2018  Bekaert David 
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
% By David Bekaert - February 2018


if strcmpi(model,'era') || strcmpi(model,'merra') || strcmpi(model,'merra2') 
    timelist = ['0000' ; '0600' ; '1200' ; '1800' ; '0000'];
elseif strcmpi(model,'era5')
    timelist = ['0000' ; '0100';'0200';'0300';'0400' ;'0500';'0600' ;'0700';'0800';'0900';'1000' ;'1100'; '1200' ;'1300';'1400';'1500';'1600' ;'1700'; '1800' ;'1900';'2000';'2100';'2200' ;'2300'; '0000'];
elseif strcmpi(model,'narr')
    fprintf('TODO \n')
    keyboard
elseif strcmpi(model_type,'gacos')
    timelist=[];
else
    error('not a supported model')
end


model_step = 24/(size(timelist,1)-1);

