function [inc]=look2inc(la,height,lat)
%LOOK2INC look angle to incidence angle
%
%    [INCIDENCE]=LOOK2INC(LOOK_ANGLE,HEIGHT,LATITUDE)
%
%    LOOK_ANGLE = Look angle (radians) (can be vector or matrix)
%    HEIGHT     = Height of satellite (m)
%    LATITUDE   = mean latitude of ground (degrees)
%
%     Copyright (C) 2015  Bekaert David - University of Leeds
%     Email: eedpsb@leeds.ac.uk or davidbekaert.com
%     With permission from Andy Hooper
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
% by Andrew Hooper, 2010

if nargin<3
   lat=40
end

lat1=lat*pi/180;

WGS84_A=6378137.0; % semimajor axis wgs84
WGS84_B=6356752.314; % semiminor axis wgs84
Re=WGS84_A*WGS84_B./sqrt(WGS84_A^2*sin(lat1).^2+WGS84_B^2*cos(lat1).^2);

a=Re+height; % Earth centre to satellite
inc=la;

[R] = fminsearch(@(p) (a^2+p^2-2*a*p*cos(mean(la(:)))-Re^2)^2,[600])
for i=1:length(la(:))
    R = fminsearch(@(p) (a^2+p^2-2*a*p*cos(la(i))-Re^2)^2,R);
    inc(i)=pi-acos((Re^2+R^2-a^2)/2/Re/R);
end

