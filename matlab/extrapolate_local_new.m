function [xy_local_box,z_box] = extrapolate_local_new(xy_local,z,x_res,y_res,n_pixels_boundary,R)
% [xy_local_box,z_box] = extrapolate_local(xy_local,z,x_res,y_res)
% Function which computes a bounding box around the data by projection of
% the circumspherece to the edges
% input:
% x_res       		Output resolution in the X-direction given in m
% y_res       		Output resolution in the Y-direction given in m
% xy_local  		A local xy grid (preferably rotated to reduce interpolation 
%                   effects for points outside the convex hull). Needs to be specified
%                   as a 2 column matrix in km.
% z                 The data observations that need to be interpolated, with each
%                   dataset repressented by a column.
% output:
% xy_local_box      The same data with in adition a bounding box 
% z_box             The same data with in additon a bounding box
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
% modifications:
% 03/12/2013    DB  Add extra option to specify the boundary extend   
% 04/12/2013    DB  RBut the box value to the closest non-NaN value

if nargin<5 || isempty(n_pixels_boundary)
    n_pixels_boundary = 5;
end
if nargin<6 || isempty(R)
    R = 0.2;        % Points within the minimum distance and this range [km] are 
                    % averaged to find the datavalue at the grid point location.
end
flag_control_fig = 0;       % when 1 plot the control figures


% converting the resolutions to km
x_res = x_res./1000;	% [km]
y_res = y_res./1000;	% [km]

% making sure the inputted data is given as a double
xy_local = double(xy_local);
z = double(z);

% computed grid extremes
[xy_min] = min(xy_local);
[xy_max] = max(xy_local);

% fixing the grid for points that will be outside the convex hull.
% An average of a region cell along the circumspherence 
x_res_area = n_pixels_boundary*x_res;
y_res_area = n_pixels_boundary*y_res;

% the vector of the bins
% Make the slightly bigger than the orgiginal grid

x_area_bins = [xy_min(1)-x_res_area:x_res_area:ceil(xy_max(1)./x_res_area).*x_res_area+x_res_area]';
y_area_bins = [xy_min(2)-y_res_area:y_res_area:ceil(xy_max(2)./y_res_area).*y_res_area+y_res_area]';


grid_bottom = [x_area_bins repmat(y_area_bins(1),length(x_area_bins),1)];
grid_left = [repmat(x_area_bins(1),length(y_area_bins)-2,1) y_area_bins(2:end-1)];
grid_top = [x_area_bins repmat(y_area_bins(end),length(x_area_bins),1)];
grid_right = [repmat(x_area_bins(end),length(y_area_bins)-2,1) y_area_bins(2:end-1)];
xy_grid = [grid_bottom ; grid_left ; grid_top ; grid_right];
clear grid_bottom grid_left grid_top grid_right


% Initializing the grid edges and their corresponding data values
n_grid_poits = size(xy_grid,1);
n_datasets = size(z,2);
z_grid = NaN([n_grid_poits n_datasets]);

% interferograms which contain only NaN vlaues
ix_ifg = find(sum(isnan(z))==n_grid_poits);
if ~isempty(ix_ifg)
   fprintf(['This is not allowed latter on in the code, as all NaNs are removed for all pixels [DB]\n'])
   keyboard
end

% removing NaN values
[row_ix,column_ix] = find(isnan(z)==1);
row_ix = unique(row_ix);
xy_local_temp =xy_local;
z_temp = z;
xy_local_temp(row_ix,:)=[];
z_temp(row_ix,:)=[];

if R==inf
    % take the average of these points to set as new grid value
    z_grid = repmat(nanmean(z_temp,1),n_grid_poits,1);   
else

    for k=1:n_grid_poits
        % compute the distance from the grid point to the dataset points
        Distance = sqrt((xy_local_temp(:,1)-xy_grid(k,1)).^2+(xy_local_temp(:,2)-xy_grid(k,2)).^2);

        % The distance to the closest point
        min_distance = min(Distance);

        % find the points within the minimum distance and a range R
        ix = find(Distance-(min_distance+R)<=0);    
        clear Distance

        % take the average of these points to set as new grid value
        z_grid(k,:) = nanmean(z_temp(ix,:),1);   



        if flag_control_fig==1 
            if k==1
                h1 = figure('name','Extrapolated data');
            end
            hold off
            plot(xy_grid(k,1),xy_grid(k,2),'r*')
            hold on
            plot(xy_local_temp([1:500:end],1),xy_local_temp([1:500:end],2),'k.')
            hold on
            scatter3(xy_local_temp([ix],1),xy_local_temp([ix],2),z_temp([ix],1),20,z_temp([ix],1),'filled')
            hold on
            view(0,90)
            axis equal
            axis tight
            pause(0.1)
        end

        clear ix

    end
end

% 
% % remove duplication when there:
% [xy_new ix temp] = unique(xy_new,'rows');
% clear temp
% z_new = z_new(ix,:);
% clear ix


% the output data
xy_local_box = [xy_local;xy_grid];
z_box = [z; z_grid];


if flag_control_fig==1 
    figure('name','Data with bounding box around it')
    scatter3(xy_local_box([1:500:length(z) length(z):end],1),xy_local_box([1:500:length(z) length(z):end],2),z_box([1:500:length(z) length(z):end],1),20,z_box([1:500:length(z) length(z):end],1),'filled')
    view(0,90)
    axis equal
    axis tight
end

%%
clear xy_local xy_grid z_grid z