function [X_regular,Y_regular,Z_regular] = interpolate_regular(xy_local,z,x_res,y_res,int_method) 
% [xy_regular_local,z_regular] = interpolate_regular(xy_local,z,x_res,y_res,int_method)
% Interpolates the data to a regular grid and fills the gaps outside the convexhull.
% The data is outputed by the new regular grid in X, in Y. The data is a regular grid Z
% with in its thrid dimentsion the the different datasets.
% input:
% x_res       		Output resolution in the X-direction given in m
% y_res       		Output resolution in the Y-direction given in m
% int_method  		The interpoaltion method to be used. This is an optional input 
%                   argument. By default a triangular interpolation is performed.
% xy_local  		A local xy grid (preferably rotated to reduce interpolation 
%                   effects for points outside the convex hull). Needs to be specified
%                   as a 2 column matrix in km.
% z                 The data observations that need to be interpolated, with each
%                   dataset repressented by a column.
% output:
% X_regular         Regular X-grid in km
% Y_regular         Regular Y_grid in km
% Z_regular         Regular Z-grid (multi-dimensional when having more than 1 dataset)
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
% modifications
% 04/12/2013    DB  do not use NaN values in the input grid during interpolation.

% setting the default interpolation method when needed
if nargin < 5 
	int_method = [];
end
if isempty(int_method)==1
	int_method = 'linear';
end
flag_control_fig = 0;       % when 1 plot the control figures


% converting the resolutions to km
x_res = x_res./1000;	% [km]
y_res = y_res./1000;	% [km]

% The number of dataset specified
n_datasets = size(z,2);

% computed grid extremes
[xy_min] = min(xy_local);
[xy_max] = max(xy_local);

% constructing a regular grid based on the resolution given 
% set-up the regualr grid
[X_regular,Y_regular] = meshgrid([xy_min(1):x_res:xy_max(1)],[xy_min(2):y_res:xy_max(2)]);
% interpolate the grid:
Z_regular = NaN([size(X_regular) n_datasets]);
tic
for k=1:n_datasets
    ix_nan = isnan(z(:,k));
	Z_regular(:,:,k) = griddata(xy_local(~ix_nan,1),xy_local(~ix_nan,2),z(~ix_nan,k),X_regular,Y_regular,int_method);
end
toc

% % Compare with interpolation of delaungy in terms of speed and results.
% Z_regular2  = NaN([size(X_regular) n_datasets]);
% tic
% DT = DelaunayTri(xy_local(:,1),xy_local(:,2));
% keyboard
% for k=1:n_datasets
%     F = TriScatteredInterp(DT,z(:,k));
%     Z_regular2(:,:,k) = F(X_regular,Y_regular);
% end
% toc


if flag_control_fig==1
    figure
    imagesc([X_regular(1,1) X_regular(1,end)],[Y_regular(1,1) Y_regular(end,1)],Z_regular(:,:,1));
    axis xy
    axis equal
    axis tight
    colorbar
end
