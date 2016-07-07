function [xy_rotatelocal,xy_local,origin] = ll2rotatelocal(ll,heading,origin)
% [xy_rotatelocal] = ll2rotatelocal(llh,heading,origin)
% Convert the longitude latitude data to a local reference frame with
% an origin in the lower left corner of the image. The heading of the satellite, 
% which is define positive clockwize from the north to the flight direction.
% The dataset will be rotated in opposit direction 
% The output grid is given in km.
%
% inputs: 
% llh:              A matrix with in its columns longitude [deg] and latitude 
%                   [deg].
% rotation_angle:   Rotation angle of the dataset [deg].
%
% Optional input:   
% origin:           lon lat of the origin [deg]
% 
% outputs:
% xy:               The data converted to a local rotated reference frame
%                   with the origin at the lower left corner.
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
% By David Bekaert - February 2013
% PhD student - University of Leeds
%
% modifications:
% 

plot_flag = 0;

if nargin<2 && nargin==1
   error('myApp:argChk', ['Heading or rotation angle needs to be specified. \n'])       
end
if nargin<1
   error('myApp:argChk', ['LonLat needs to be specified. \n'])       
end
if nargin<3
    % set the origin of the dataset
    origin = mean(ll(:,[1 2]));                % [deg]
end
    

% convert from geo-coordiantes to local reference frame
xy_local = llh2local(ll(:,[1 2])',origin);

% define and convert the rotation angle to radians:
rotation_angle = -heading*pi/180;     % [rad]

% Rotating the local reference frame
Rotation_matrix = [cos(rotation_angle)  -sin(rotation_angle)
                   sin(rotation_angle)  cos(rotation_angle)];
% performing the rotation of the data
xy_rotatelocal = (Rotation_matrix\xy_local)';     

% putting the local data in right format
xy_local = xy_local';


if plot_flag==1 
    figure('name','Dataset in geo-coordinates')
    plot(ll(:,1),ll(:,2),'k.')
    axis equal
    axis tight

    figure('name','Rotated dataset in a local reference frame')
    plot(xy_rotatelocal(:,1),xy_rotatelocal(:,2),'k.')
    axis equal
    axis tight
end


