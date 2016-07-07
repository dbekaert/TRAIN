function [z] = griddata_version_control(gridX,gridY,gridZ,x,y,method,matlab_version)
% This function just computes the gridded data set, but its input is
% modified to work with matlab version newer than 2012.
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
% By David Bekaert - December 2012

if matlab_version>=2012
    z=griddata(gridX,gridY,gridZ,x,y,method);
else
    z=griddata(gridX,gridY,gridZ,x,y,method,{'QJ'});
end
