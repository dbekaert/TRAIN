function [z_local] = interpolate2local(X_regular,Y_regular,Z_regular,xy_local,int_method) 
% [z_local] = interpolate2local(X_regular,Y_regular,Z_regular,xy_local,int_method) 
% Interpolates thed ata to a refular grid and fill the gaps outside the convexhull.
% The data is outputed by the new regular grid in X, in Y. The data is a regular grid Z
% with in its thrid dimentsion the the different datasets.
% input:
% X_regular         Regular X-grid in km
% Y_regular         Regular Y_grid in km
% Z_regular         Regular Z-grid (multi-dimensional when having more than 1 dataset)
% int_method  		The interpoaltion method to be used. This is an optional input 
%                   argument. By default a triangular interpolation is performed.
% xy_local  		A local xy grid (preferably rotated to reduce interpolation 
%                   effects for points outside the convex hull). Needs to be specified
%                   as a 2 column matrix in km.
% output:
% z_local           The data observations that have been interpolated, with each
%                   dataset repressented by a column.
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

% setting the default interpolation method when needed
flag_control_fig = 0;       % when 1 plot the control figures
if nargin < 5 
	int_method = [];
end
if isempty(int_method)==1
	int_method = 'linear';
end

% Convert the position grid back to a column vector
x_regular = reshape(X_regular,[],1);
clear X_regular
y_regular = reshape(Y_regular,[],1);
clear Y_regular

% Converting the datagrid back to a column vector or a matrix 
n_datasets = size(Z_regular,3);
n_gridpoints = length(y_regular);
z_regular = reshape(Z_regular,n_gridpoints,n_datasets,1);

z_local = NaN([size(xy_local,1) n_datasets]);
for k=1:n_datasets
%     % alternative replace by tri interpolation
%     F = TriScatteredInterp(x_regular,y_regular,z_regular(:,k));
%     z_local_tri(:,k) = F(xy_local);
	z_local(:,k) = griddata(x_regular,y_regular,z_regular(:,k),xy_local(:,1),xy_local(:,2),int_method);
end

if flag_control_fig==1
    figure
    scatter3(xy_local(:,1),xy_local(:,2),z_local(:,1),3,z_local(:,1),'filled');
    hold on
    scatter3(x_regular(:,1),y_regular(:,1),z_regular(:,1),50,z_regular(:,1),'filled');
    view(0,90)
    axis xy
    axis equal
    axis tight
    colorbar
    clear x_regular y_regular z_regular
end