function [xyz_input,xyz_output] = load_wrf_SAR(filename,xy_out_grid,poly_crop_filename) 
% [xyz] = load_wrf_SAR(filename,crop_filename,geocoord_flag) 
% script to read and optional interpolate or crop the data. When
% interpolating, an initial crop is done based on the extend of the to be
% interpolated coordinated. It is assumed that locations are
% given in geo-coordinates (degrees).
%
% INPUTS:
% filename              Filename with full path to the ".xyz" dataset to be read
% 
% OPTIONAL INPUTS
% xy_out_grid           Grid to which the input data should be interpolated
%                       to. This should have the same coordinates as the
%                       input dataset. A triangular interpolation is
%                       assumed.
% poly_crop_filename    Filename to the .mat crop file. The coordinates needs
%                       be the same as the xy coordinates of the datafile. The
%                       variable in the poly_crop_filename is needs to be called 
%                       "poly.xy".
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
% By David Bekaert - University of Leeds

% This is the same for ERA-I
if nargin <3 || isempty(poly_crop_filename)
   poly_crop_filename =[];
end
if nargin <2 || isempty(xy_out_grid)
   xy_out_grid =[];
end

[xyz_input,xyz_output] = load_era_SAR(filename,xy_out_grid,poly_crop_filename);
