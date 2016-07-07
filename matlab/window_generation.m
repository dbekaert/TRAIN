function [iterate] = window_generation(xy_local,n_patches)
% 
%  Windows are increasing bottom up and start in the lower left corner.
%
% Inputs:
% xy_local      Local grid, recomended to be rotated to reduce number of windows
%               outside the dataset, specified as a 2 column matrix in km.
% Optional inputs:
% iterate           A struct containign vairables that are computed for the
%                   first run and whicha re not recompeted for the other runs.
%                   This includes:
%                   window_ix: a variable struct containing which points 
%                               are in wich window.
%                   window_xy: the window positions that cover the data
%                   window_xy_extra: the window positions of the bounding box 
%                               around the data.
% crop_flag     	Set to 'y' to crop out an area. By default this is not done.
%               	When 'y': Filename = "area_ex.mat", with variable called lonlat.
%               	Use: lonlat=ginput; to select polygon around the to be cropped (out) area.
%               	Save as: save('area_ex.mat','lonlat');
%
%     window_generation.m
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
%
% September 2010 --- Bekaert David --- Initial codings
%
% modifications:
% 07/2014   DB:     Fix the cropping out of a region

%% defaults
% Defining input variables   
if nargin<1
    error('myApp:argChk', ['Too few input arguments...\nAbort... \n'])    
end

 


flag_figure_checks =0;
n_points_min = 1;              % minimum number of points in a cell otherwize rejected

patch_overlap = 50;




% getting the parameters from the parm list
 
crop_flag = getparm_aps('crop_flag'); 
stamps_processed = getparm_aps('stamps_processed');



%% Redefine the grid wrt lower left corner.
xy_local(:,2) = xy_local(:,2) - min(xy_local(:,2));
xy_local(:,1) = xy_local(:,1) - min(xy_local(:,1));


%% Computing the patch edges
% only perform operation this when it has not been done before
% size of the image in km
width_x = (max(xy_local(:,1))-min(xy_local(:,1)));
width_y = (max(xy_local(:,2))-min(xy_local(:,2)));  

% Assuming the patches to be approximately square to find number of 
% patches in x and y direction.
patch_width = sqrt(width_x*width_y./n_patches);
n_patches_x = ceil(width_x./patch_width);
n_patches_y = ceil(width_y./patch_width);
fprintf([num2str(n_patches_x),' patches (+2 for the edges) in x-dir \n']);
fprintf([num2str(n_patches_y),' patches (+2 for the edges) in y-dir \n']);

% Assume same patch overlap
% But allow in futuer for a different one in x and y direction.
patch_overlap_x = patch_overlap;
patch_overlap_y = patch_overlap;
fprintf([num2str(patch_overlap_x),' procent patch overlap in x direction \n']);
fprintf([num2str(patch_overlap_y),' procent patch overlap in y direction \n']);

% size of the patches including overlap
patch_width_x = width_x./(patch_overlap_x./100 + (1-patch_overlap_x./100)*n_patches_x);
patch_width_y = width_y./(patch_overlap_y./100 + (1-patch_overlap_y./100)*n_patches_y);

fprintf(['Patch width with overlap in x-direction: ' , num2str(patch_width_x),'\n']);
fprintf(['Patch width with overlap in y-direction: ' ,num2str(patch_width_y),'\n']);

% Patch size center to center
patch_c2c_y = patch_width_y - patch_width_y*patch_overlap_y./100;
patch_c2c_x = patch_width_x - patch_width_x*patch_overlap_x./100;
patch_c2c = max([patch_c2c_y patch_c2c_x]);                         % size of the gaussian filter for weigthing on distance [km]        



%% determining the points in each window 
% Paches are increasing bottom up and start in the lower left corner.
% only perform the following operation when it has not been done before.
% doing onyl onces will save time for other interferograms
n_points = size(xy_local,1);
ix_window = [];             % number of the selected windows with respect to the intial set of generated windows.
ixy_local = [[1:n_points]' xy_local];
counter = 1;                    % counter of the windows with points
counter_all = 1;                % counter of all windows including those without points
for i=1:n_patches_x+2
    % search for the points in the window
    x_bounds    = [(i-1)*patch_width_x*(1-patch_overlap_x./100) (i-1)*patch_width_x*(1-patch_overlap_x./100)+patch_width_x]-patch_width_x/2;
    ix_temp = find(x_bounds(1)<= ixy_local(:,2) &  ixy_local(:,2)< x_bounds(2));
    ixy_local_temp = ixy_local(ix_temp,:);
    clear ix_temp
    for k=1:n_patches_y+2
        % search for the points in the window
        y_bounds  = [(k-1)*patch_width_y*(1-patch_overlap_y./100) (k-1)*patch_width_y*(1-patch_overlap_y./100)+patch_width_y]-patch_width_y/2;
        ix_temp = find(y_bounds(1)<= ixy_local_temp(:,3) & ixy_local_temp(:,3)< y_bounds(2));
        ix = ixy_local_temp(ix_temp,1);
        clear ix_temp 

        % keeping only those windows were there is actual data, and
        % were there are more than n_points_min
        if isempty(ix)~=1 && length(ix)>=n_points_min
            % saving the data back into the struct
            window_ix{counter} = ix;

            % computing the center of the window
            window_xy(counter,1) = x_bounds(1)+(x_bounds(2)-x_bounds(1))/2;
            window_xy(counter,2) = y_bounds(1)+(y_bounds(2)-y_bounds(1))/2;

            % plotting the different patches
            if flag_figure_checks==1
                patch_colors = jet((n_patches_x+2)*(n_patches_y+2));
                if counter==1
                    h1 = figure('name','grid');
                else
                    figure(h1)
                    hold on
                end
                plot(window_xy(counter,1),window_xy(counter,2),'o','color',patch_colors(counter,:))
                hold on
                plot([x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1) x_bounds(1)],[ y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1)],'-','color',patch_colors(counter,:))
                axis equal 
                axis tight
                xlim([0 width_x])
                ylim([0 width_y])
            end


            % saving the orignal window number for the windows with points
            ix_window(counter)=counter_all;
            counter = counter+1;

        end
        counter_all = counter_all+1;
        clear ix
    end
end
clear x_bounds y_bounds counter


% the windows that are edges:
edges = [1:n_patches_y+2 (n_patches_y+2)*(n_patches_x+2)-(n_patches_y+2)+1:(n_patches_x+2)*(n_patches_y+2) 1:n_patches_y+2:(n_patches_y+2)*(n_patches_x+2)-n_patches_y+1   n_patches_y+2:n_patches_y+2:(n_patches_y+2)*(n_patches_x+2)];
edges = unique(edges);
% search those windows with points that are edge windows
[edges,temp,ix_edges_window] = intersect(edges,ix_window);
clear counter_all edges temp


% getting the final number of windows
n_windows = length(window_ix);

% Computing weights based on distance from the points to the windows centers
X_points = repmat( xy_local(:,1),1,n_windows);
Y_points = repmat(xy_local(:,2),1,n_windows);
X_windows = repmat(window_xy(:,1)',n_points,1);
Y_windows = repmat(window_xy(:,2)',n_points,1);
Distance = sqrt((X_points-X_windows).^2+(Y_points-Y_windows).^2);
clear X_points Y_points X_windows Y_windows
% Compute the gaussian weight based on distance

% 
% w_d = normpdf(Distance,0,patch_c2c);                               % [n_points x n_windows]
% if flag_figure_checks==1
%     figure('Name','Weighting based on Distance')
%     Distance_temp = [min(min(Distance)):(max(max(Distance))-min(min(Distance)))/100:max(max(Distance))]';
%     weights_temp = normpdf(Distance_temp,0,patch_c2c); 
%     plot(Distance_temp,weights_temp,'b.-')
% end
% clear Distance  


% getting the variables in a structure to load on the next run
iterate.window_ix = window_ix;
iterate.window_xy = window_xy;
iterate.ix_edges_window = ix_edges_window;
% iterate.w_d = w_d;



