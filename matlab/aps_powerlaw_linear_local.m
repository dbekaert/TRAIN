function [ph_tropo_powerlaw,ph_tropo_powerlaw_band,iterate,slope_local,slope_local_band,perc_cons_K] = aps_powerlaw_linear_local(xy_local,local_origin,data_local1,data_local2,iterate,dimension_filter,save_path_outlier_window,ifg_number,ix_phase_nan_interferogram,debug_fig_linear)
% 
% Function that computes the linear relation between data_local1 and
% data_local2, assuming the slope is computed with data_local1 on the x-axis 
% data_local2 on the y-axis. Windows are increasing bottom up and start in
% the lower left corner.
%
% Inputs:
% xy_local      Local grid, recomended to be rotated to reduce number of windows
%               outside the dataset, specified as a 2 column matrix in km.
% data_local1   Dataset 1 (matrix with n column bandfiltered datasets) X-axis
% data_local2   Dataset 2 (matrix with n column bandfiltered datasets) Y-axis
% 
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
%
% OPTIONAL: visulize the estimation for an interferogram and a specific
% patch by changing the default variable "test" to 1 and "test_patch" to the
% specific patch number (counted columswize from top left in radar coordinates). 
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
%
% September 2010 --- Bekaert David --- Initial codings
%
% 12/2012:  DB    Start from previous version of the code
% 01/2013:  DB    Allow for a wider range of bandfiltered data to be inputed.
% 02/2013:  DB    Adding estimation of slope from multiple bands
% 02/2013:  DB    Fix a bug when running for multiple datasets, i.e. when already 
%                 part of the window info is feeded directly into the function.
% 03/2013:  DB	  Allow for an extra argument (filepath), which when specified is
%                 used to save the results of the outlier rejection for each window.
% 03/2013:  DB	  Specify the number of patches including those for which the estimation
%                 is based on half a patch. These were introduced to reduce edge effects.
% 03/2013:  DB    Allow for an extra input argument that specifies the type
%                 of band filtering that was applied and displays it in the figures.
% 04/2013   DB:   Estimation of the slope at local grid by gaussian kernel
% 04/2013   DB:   Removing the extra input parameters and load them from the parm file
% 04/2013   DB:   Allow for the computation of the delay for each bandfilter
% 10/2013   DB:   Adding the option to have a single patch (n_patch=0)
%                 and using a planar estiamte for K
% 12/2013   DB:   Give the local slope also as output 
% 02/2014   DB:   Output the local slope for the bands too.
% 02/2014   DB:   Include a flag that allows for ridge-topo information to
%                 be used to constrain patching
% 05/2014   DB:   Fix dummy output when there is none.
% 07/2014   DB:   Include a cropped region for the patches option. Has been
%                 tested for ridges and when more than 1 window is left
% 06/2015   DB:   Change the default of the planar flag to 'n'
% 01/2016   DB:   Include the planar mode as parameter, fix for non-stamps
%                 processed data
% 02/2016   DB:   Include a variable which track NaN pixels and which
%                 updates the variable of slope_dataset

%% defaults
% Defining input variables   
if nargin<4
    error('myApp:argChk', ['Too few input arguments...\nAbort... \n'])    
end

flag_figure_checks =0;
debug_fig = 0;
n_points_min = 50;      % minimum number of points in a cell otherwize rejected
n_boot_runs =1000;     % number of bootstrap runs uised to estimate the uncertaincy
bootstrap_flag = 1;     % when 1 uses bootstrapping for the uncertaincy estimation. otherwize its propagation of errors.


% Get matlab version as function arguments change with the matlab version
matlab_version = version('-release');           % [DB] getting the matlab version
matlab_version = str2num(matlab_version(1:4));  % [DB] the year


if nargin <5 || isempty(iterate)==1
    window_ix=[];  % force the window to be re-estimated
    window_xy=[];
    window_xy_extra = [];
    ix_edges_window = [];
    compute_windows = 1;
else
    compute_windows = 0;
end
if nargin <6 || isempty(dimension_filter)
    dimension_filter = [];
end
if isempty(dimension_filter)~=1 && length(dimension_filter)~=size(data_local2,2)
   fprintf(['Warning: dimension_filter has different size then the dataset and will be ignored during plotting. \n']) 
   dimension_filter=[];
end

if nargin<7 || isempty(save_path_outlier_window)
    save_flag_outlier_windows = 0;
else
    save_flag_outlier_windows = 1;
end
if nargin<8
    fprintf(['Assuming the power-law has been applied using one alpha and h0 \n']);
    ifg_number=1;
end
if nargin<9
    ix_phase_nan_interferogram = [];
end

if nargin<10
    debug_fig_linear=0;
end

fontsize = 15;


% getting the parameters from the parm list
n_patches = getparm_aps('powerlaw_n_patches');
powerlaw_all_bands = getparm_aps('powerlaw_all_bands');
powerlaw_ridge_constraint = getparm_aps('powerlaw_ridge_constraint');
patch_overlap = getparm_aps('powerlaw_patch_overlap'); 
crop_flag = getparm_aps('crop_flag'); 
h0 = getparm_aps('powerlaw_h0');
alpha = getparm_aps('powerlaw_alpha');
stamps_processed = getparm_aps('stamps_processed');
plane_mode = getparm_aps('powerlaw_plane_mode');       % run in a planar mode correction in case the whole study area is contained in one window

if strcmp(stamps_processed,'y')
    load psver
else
    psver = 2; 
end
 % loading the topography information
hgt_matfile = getparm_aps('hgt_matfile');
hgt = load(hgt_matfile);
hgt = hgt.hgt; 


%% Redefine the grid wrt lower left corner.
xy_local_regrid = [min(xy_local(:,1)) min(xy_local(:,2))];
xy_local(:,2) = xy_local(:,2) - xy_local_regrid(:,2);
xy_local(:,1) = xy_local(:,1) - xy_local_regrid(:,1);



%% Computing the patch edges
% only perform operation this when it has not been done before
if compute_windows==1  && strcmpi(powerlaw_ridge_constraint,'n')
    
    
    % size of the image in km
    width_x = (max(xy_local(:,1))-min(xy_local(:,1)));
    width_y = (max(xy_local(:,2))-min(xy_local(:,2)));  

    % getting the variables from the parm_aps file
    heading = getparm_aps('heading');

    if n_patches~=0
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
    else
       fprintf('Use an individual window \n') 
       patch_width_x=width_x;
       patch_width_y=width_y;
       fprintf(['Patch width with overlap in x-direction: ' , num2str(width_x),'\n']);
       fprintf(['Patch width with overlap in y-direction: ' ,num2str(width_y),'\n']);
    end
    
    powerlaw_window_size = [patch_width_x patch_width_y];
    if exist('tca_support.mat')==2
       save(['tca_support.mat'],'-append','powerlaw_window_size') 
    else
       save(['tca_support.mat'],'powerlaw_window_size') 
    end
    
    
elseif compute_windows==1  && strcmpi(powerlaw_ridge_constraint,'y')  
    tca_support = load('tca_support.mat','powerlaw_ridges');
    tca_support = tca_support.powerlaw_ridges;
    mountain_ridge = tca_support.mountain_ridge;
    hard_ridge_flag = tca_support.hard_ridge_flag;
    
    mountain_ridge =  mountain_ridge(~cellfun('isempty',mountain_ridge));
    mountain_ridge{length(mountain_ridge)+1}=tca_support.InSAR_convexhull;
    hard_ridge_flag = [hard_ridge_flag 0];
    
    % getting the variables from the parm_aps file
    heading = getparm_aps('heading');
    
    % loading the lonlat information 
    xy_rotatelocal_all=[];
    min_region=[];
    max_region=[];
    for kk=1:length(mountain_ridge)
        % rotating the data
        [xy_rotatelocal{kk},temp{kk}] = ll2rotatelocal(mountain_ridge{kk},heading,local_origin);
        xy_rotatelocal{kk}=round(xy_rotatelocal{kk});
        max_region = max([xy_rotatelocal{kk} ; max_region]);
        min_region = min([xy_rotatelocal{kk} ; min_region]);
        clear temp
        
        % storing the InSAR bounding box for latter use
        if kk==length(mountain_ridge)
            xy_InSAR_convexhull = xy_rotatelocal{kk};
        end
    end
     
    
    matrix_patches = zeros([max_region - min_region]+1);
    y_vector=[];
    x_vector = [];
    z_vector = [];
    for kk=1:length(mountain_ridge)
        
        if hard_ridge_flag(kk)==1
            z_ridge = 10000;
        else
            z_ridge = 100;
        end
        
        clear ix_lines  
        % shift the origin of the indexes
        ix_lines(:,1) = xy_rotatelocal{kk}(:,1) - min_region(1,1)+1;
        ix_lines(:,2) = xy_rotatelocal{kk}(:,2) - min_region(1,2)+1;
    
        for kkk=1:length(ix_lines(:,1))-1 
            x_sampling =  [min([ix_lines(kkk,1) ix_lines(kkk+1,1)]): max([ix_lines(kkk,1) ix_lines(kkk+1,1)])]';
            if length(x_sampling)<=1
                % The x points are at the same coordinate, i.e. this is a
                % vertical line, change the sampling to be in y-direction.
                y_sampling =  [min([ix_lines(kkk,2) ix_lines(kkk+1,2)]): max([ix_lines(kkk,2) ix_lines(kkk+1,2)])]';
                if length(y_sampling)>1
                     clear x_sampling
                     x_sampling = interp1(ix_lines(kkk:kkk+1,2),ix_lines(kkk:kkk+1,1),y_sampling);
                     x_vector = [x_vector ;x_sampling];
                     y_vector = [y_vector ;y_sampling];
                     z_vector = [z_vector ;z_ridge.*ones(size(y_sampling))];
                else
                    % this is a single point -- appears to occur double in
                    % the line.
                    x_vector = [x_vector ;x_sampling(1)];
                    y_vector = [y_vector ;y_sampling(1)];
                    z_vector = [z_vector ;z_ridge.*ones(size(y_sampling))];
                end

            else
                % do the interpolation in x-direction
                if ~isempty(x_sampling)
                    y_sampling = interp1(ix_lines(kkk:kkk+1,1),ix_lines(kkk:kkk+1,2),x_sampling);
 
                     x_vector = [x_vector ;x_sampling];
                     y_vector = [y_vector ;y_sampling];
                     z_vector = [z_vector ;z_ridge.*ones(size(y_sampling))];

        

                     if abs(mean(diff(y_sampling)))>1.5
                        % step in y direction is large -- change
                        % interpolation direction to scope with this
                         y_sampling = [min([ix_lines(kkk,2) ix_lines(kkk+1,2)]): max([ix_lines(kkk,2) ix_lines(kkk+1,2)])]';
                         clear x_sampling
                         x_sampling = interp1(ix_lines(kkk:kkk+1,2),ix_lines(kkk:kkk+1,1),y_sampling);
                         x_vector = [x_vector ;x_sampling];
                         y_vector = [y_vector ;y_sampling];
                         z_vector = [z_vector ;z_ridge.*ones(size(y_sampling))];
                     
                    end
                end
            end
            clear x_sampling y_sampling
        end
        clear z_ridge
    end
    x_vector = round(x_vector) ;
    y_vector = round(y_vector) ;
    
    ix = sub2ind(size(matrix_patches), x_vector, y_vector);
    matrix_patches(ix)=z_vector;
    matrix_patches = matrix_patches';

    
    % extend the lines to remove holes
    H = fspecial('average',2);
    matrix_patches_mean = imfilter(matrix_patches,H,'replicate');
    matrix_patches_mean(matrix_patches_mean>=1000)=10000;
    matrix_patches_mean(matrix_patches_mean>0 & matrix_patches_mean<1000)=100;

    clear H
    
    % shifting the grid axis to the same as that of the InSAR points
    x_temp= [min_region(1):max_region(1)]- xy_local_regrid(:,1);
    y_temp = [min_region(2):max_region(2)]- xy_local_regrid(:,2);
    % putting hte InsAR convechull in the same reference frame
    xy_InSAR_convexhull(:,1) = xy_InSAR_convexhull(:,1)- xy_local_regrid(:,1);    
    xy_InSAR_convexhull(:,2) = xy_InSAR_convexhull(:,2)- xy_local_regrid(:,2);
    
    
    % generating a mesh grid 
    [Xmesh_patch_regions,Ymesh_patch_regions] = meshgrid(x_temp,y_temp);

    % watershed computation
    patch_regions = watershed(matrix_patches_mean);
    patch_regions = double(patch_regions);
    
    % remove the watershed edge and replace it to the patch value the closest by
    [ix_row ix_column ] = find(patch_regions==0);
    ix_zeros = sub2ind(size(patch_regions), ix_row, ix_column);
    
    % coordinates for the window edges:
    xy_edges = [ix_column ix_row]+repmat(xy_local_regrid,size(ix_column,1),1);
    [window_edges_ll] = rotatelocal2ll(xy_edges,heading,local_origin);
    
    % remove the watershed edges for interpoaltion
    Xmesh_patch_regions_temp = reshape(Xmesh_patch_regions,[],1);
    Ymesh_patch_regions_temp = reshape(Ymesh_patch_regions,[],1);
    patch_regions_temp = reshape(patch_regions,[],1);
    patch_regions_temp(ix_zeros)=[];
    Xmesh_patch_regions_temp(ix_zeros) = [];
    Ymesh_patch_regions_temp(ix_zeros) = [];
    
    
    [z] = griddata_version_control(Xmesh_patch_regions_temp,Ymesh_patch_regions_temp,patch_regions_temp,Xmesh_patch_regions(ix_zeros),Ymesh_patch_regions(ix_zeros),'nearest',matlab_version);
    patch_regions(ix_zeros)=z;
    clear z
        

    % removing those watershed that are not part of the insar region
    patch_regions_temp = reshape(patch_regions,[],1);
    Xmesh_patch_regions_vector = reshape(Xmesh_patch_regions,[],1);
    Ymesh_patch_regions_vector = reshape(Ymesh_patch_regions,[],1);
    IN = inpolygon(Xmesh_patch_regions_vector,Ymesh_patch_regions_vector,xy_InSAR_convexhull(:,1),xy_InSAR_convexhull(:,2));
    patch_regions_temp(IN)=[];

    % potential patches outside the insar region
    non_patches_pot = unique(patch_regions_temp);
    non_patches_pot(isnan(non_patches_pot))=[];
    
    
    % check if there are enough counts in those windows outside the insar
    % region, as there might eb a slight grid offset. Otherwize keep the
    % patch as a valid patch.
    non_patches = [];
    for kk=1:length(non_patches_pot)
        count = sum(non_patches_pot(kk)==patch_regions_temp);
        count_total = sum(sum(non_patches_pot(kk)==patch_regions));
        if (count./count_total)*100>10
            non_patches = [non_patches ; non_patches_pot(kk)];    
        end
    end
    % removing those bad located patched from the watershed
    for kk=1:length(non_patches)
        patch_regions(patch_regions==non_patches(kk))=NaN;
    end

    % redefine the counting of the patches from 1 to number of patches
    patch_id = unique(patch_regions);
    patch_id(isnan(patch_id))=[];

    n_patches = length(patch_id);
    for k=1:n_patches
       patch_regions(patch_regions==patch_id(k))=-2*n_patches+k;
    end
    patch_regions=patch_regions+2*n_patches;
    
    
    if flag_figure_checks==1 && debug_fig==1

        figure
        plot(xy_local(:,1),xy_local(:,2),'k.')
        axis equal
        axis tight
        hold on
        plot(xy_InSAR_convexhull(:,1),xy_InSAR_convexhull(:,2),'r-')

        figure
        imagesc(matrix_patches)
        axis equal
        axis tight
        axis xy

        figure
        imagesc(matrix_patches_mean)
        axis equal
        axis tight
        axis xy

        figure
        imagesc(patch_regions)
        axis equal
        axis tight
        axis xy
    end   
    
        
end


%% determining the points in each window 
% only perform patch point computation when it has not been done before.
% doing only onces will save time for other interferograms
n_points = size(data_local1,1);
ix_window = [];             % number of the selected windows with respect to the intial set of generated windows.
if compute_windows==1   && strcmpi(powerlaw_ridge_constraint,'n')
    % this is the patch approach.
    % Paches are increasing bottom up and start in the lower left corner.
    if n_patches==0
        % take only an individual window
        window_ix{1}=[1:n_points]';
        window_xy(1,1)  = min(xy_local(:,1)) + width_x/2;
        window_xy(1,2)  = min(xy_local(:,2)) + width_y/2;
        window_xy_extra = [];
        ix_edges_window = [];
        w_d = ones([n_points 1]);
    else
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

                    % keeping the information of the window location such
                    % it can be converted to lon lat later on;
                    window_box{counter} = [[x_bounds(1);x_bounds(2);x_bounds(2);x_bounds(1);x_bounds(1)],[ y_bounds(1);y_bounds(1);y_bounds(2);y_bounds(2);y_bounds(1)]];
                    
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

        % convert the windows boxes into lon lats and save them for support
        % information
        for kkk=1:length(window_box)
            window_box_temp =  window_box{kkk}+repmat(xy_local_regrid,size(window_box{kkk},1),1);
            [window_box_ll{kkk}] = rotatelocal2ll(window_box_temp,heading,local_origin);
        end
        window_xy_temp = window_xy+repmat(xy_local_regrid,size(window_xy,1),1);
        [window_box_center_ll] = rotatelocal2ll(window_xy_temp,heading,local_origin);

        
        % saving in a struct to separate from other techniques
        powerlaw_windows.window_box_ll  = window_box_ll;
        powerlaw_windows.window_box_center_ll = window_box_center_ll;
                
%         keyboard
        if exist('tca_support.mat')==2
           save(['tca_support.mat'],'-append','powerlaw_windows') 
        else
           save(['tca_support.mat'],'powerlaw_windows') 
        end
        clear window_box_ll

        % the windows that are edges:
        edges = [1:n_patches_y+2 (n_patches_y+2)*(n_patches_x+2)-(n_patches_y+2)+1:(n_patches_x+2)*(n_patches_y+2) 1:n_patches_y+2:(n_patches_y+2)*(n_patches_x+2)-n_patches_y+1   n_patches_y+2:n_patches_y+2:(n_patches_y+2)*(n_patches_x+2)];
        edges = unique(edges);
        % search those windows with points that are edge windows
        [edges,temp,ix_edges_window] = intersect(edges,ix_window);
        clear counter_all edges temp

        powerlaw_windows_edge = ix_edges_window; 
        if exist('tca_support.mat')==2
           save(['tca_support.mat'],'-append','powerlaw_windows_edge') 
        else
           save(['tca_support.mat'],'powerlaw_windows_edge') 
        end
        clear window_box_ll

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

        
        w_d = normpdf(Distance,0,patch_c2c);                               % [n_points x n_windows]
        if flag_figure_checks==1
            figure('Name','Weighting based on Distance')
            Distance_temp = [min(min(Distance)):(max(max(Distance))-min(min(Distance)))/100:max(max(Distance))]';
            weights_temp = normpdf(Distance_temp,0,patch_c2c); 
            plot(Distance_temp,weights_temp,'b.-')
        end
        clear Distance  
    end

    % getting the variables in a structure to load on the next run
    iterate.window_ix = window_ix;
    iterate.window_xy = window_xy;
    iterate.window_xy_extra = window_xy_extra;
    iterate.ix_edges_window = ix_edges_window;
    iterate.w_d = w_d;
    % next are the values in case of ridge approach. Put them to empty so
    % they will be ignored on the next run
    ridges_flag = 'n';
    patch_regions = [];
    Ymesh_patch_regions = [];
    Xmesh_patch_regions = [];
    
    iterate.ridges_flag=ridges_flag;
    iterate.patch_regions=patch_regions;
    iterate.Ymesh_patch_regions = Ymesh_patch_regions;
    iterate.Xmesh_patch_regions = Xmesh_patch_regions;


elseif compute_windows==1 && strcmpi(powerlaw_ridge_constraint,'y')
    % this is the ridges approach
    
    % in the previous section patch_regions where defined.
    % here all the PS poitns are interpolated based onthe patch region,
    % giving a patch id for point.
    z_points = interp2(Xmesh_patch_regions,Ymesh_patch_regions,patch_regions,xy_local(:,1),xy_local(:,2),'nearest');

    if flag_figure_checks==1
        figure
        scatter3(xy_local(:,1),xy_local(:,2),z_points,3,z_points,'filled')
        view(0,90)
        axis equal
        axis tight
    end
       
    % construct the iterate variable and populate the data
    for kk=1:n_patches
        ix = find(z_points==kk);
         window_ix{kk} = ix;
         window_xy(kk,:) = mean(xy_local(ix,:)); 
    end
    
    % coordinates for the window edges:
    window_xy_centers =  window_xy+repmat(xy_local_regrid,size(window_xy,1),1);
    [window_center_ll] = rotatelocal2ll(window_xy_centers,heading,local_origin);
    clear window_xy_centers
    
    if flag_figure_checks==1 
        figure('name','Patch definition') 
        plot(window_edges_ll(:,1),window_edges_ll(:,2),'k.')
        hold on
        for kk=1:size(window_center_ll,1)
            text(window_center_ll(kk,1),window_center_ll(kk,2),[num2str(kk)],'color','r','fontsize',fontsize)
            plot(window_center_ll(kk,1),window_center_ll(kk,2),'r.') 

        end
        axis equal
        axis tight
        xlabel('longitude','fontsize',fontsize)
    end
    
    
    % getting the final number of windows
    n_windows = length(window_ix);

    % comnputing the average distance between windows
    for kk=1:n_windows
        X_windows = repmat(window_xy(:,1)',n_windows,1);
        Y_windows = repmat(window_xy(:,2)',n_windows,1);
        Distance_windows = sqrt((X_windows-X_windows').^2+(Y_windows-Y_windows').^2);
        Distance_windows_sorted = sort(Distance_windows,1);
        patch_c2c = 2/3*mean(Distance_windows_sorted(2,:));
    end
    
    
    % Computing weights based on distance from the points to the windows centers
    X_points = repmat( xy_local(:,1),1,n_windows);
    Y_points = repmat(xy_local(:,2),1,n_windows);
    X_windows = repmat(window_xy(:,1)',n_points,1);
    Y_windows = repmat(window_xy(:,2)',n_points,1);
    Distance = sqrt((X_points-X_windows).^2+(Y_points-Y_windows).^2);
    clear X_points Y_points X_windows Y_windows
    
  
    dimensions = ceil(max(xy_local))-floor(min(xy_local));
    grid_resolution = 10;           % grid resolution in km;


    [Xgrid,Ygrid] = meshgrid([-1*ceil(patch_c2c./grid_resolution):grid_resolution:dimensions(1)+ceil(patch_c2c./grid_resolution)],[-1*ceil(patch_c2c./grid_resolution):grid_resolution:dimensions(2)+ceil(patch_c2c./grid_resolution)]); 
    clear dimensions
    dimensions_grid = size(Xgrid);
    Xgrid_vector = reshape(Xgrid,[],1);
    Ygrid_vector = reshape(Ygrid,[],1);


    X_points = repmat( Xgrid_vector,1,n_windows);
    Y_points = repmat(Ygrid_vector,1,n_windows);
    X_windows = repmat(window_xy(:,1)',length(Ygrid_vector),1);
    Y_windows = repmat(window_xy(:,2)',length(Ygrid_vector),1);
    Distance_grid = sqrt((X_points-X_windows).^2+(Y_points-Y_windows).^2);
    clear X_points Y_points X_windows Y_windows

    fprintf('Checking if there are hard ridges \n')
    if sum(hard_ridge_flag)>0  
%         keyboard
%         ix_hard_ridges = find(hard_ridge_flag==1);
%         xy_hard_ridge = [];
%         for kkk=1:sum(hard_ridge_flag)
%             ix_hard_ridge = ix_hard_ridges(kkk);
% 
%             % getting the coordinates of the hard ridge
%             xy_hard_ridge_temp = xy_rotatelocal{ix_hard_ridges(kkk)}-repmat(xy_local_regrid,size(xy_rotatelocal{ix_hard_ridges(kkk)},1),1);
%             xy_hard_ridge = [xy_hard_ridge ; xy_hard_ridge_temp];
%             clear xy_hard_ridge_temp
%             if kkk+1<=sum(hard_ridge_flag)
%                 xy_hard_ridge = [xy_hard_ridge; NaN NaN];
%             end
%             
%         end
%            
%         % loopign for every grid point
%         tic
%         for zz=1:length(Ygrid_vector)
% 
%             for kkkk=1:length(window_xy)
% 
%                 [x0,y0] = intersections(xy_hard_ridge(:,1),xy_hard_ridge(:,2),[window_xy(kkkk,1);Xgrid_vector(zz)],[window_xy(kkkk,2);Ygrid_vector(zz)],'robust');
%                 
%                 if ~isempty(x0)
%                     
%                     Distance_grid(zz,kkkk)=NaN;
%                 end
%                 clear xo yo
% 
%                 %INTERSECTIONS Intersections of curves.
%                 %   Computes the (x,y) locations where two curves intersect.  The curves
%                 %   can be broken with NaNs or have vertical segments.
%                 %
%                 % Example:
%                 %   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%                 %
%                 % where X1 and Y1 are equal-length vectors of at least two points and
%                 % represent curve 1.  Similarly, X2 and Y2 represent curve 2.
%                 % X0 and Y0 are column vectors containing the points at which the two
%                 % curves intersect.
%             end
%         end
%         toc

        % see if there are patches that are separated by windows.
        ix_hard_ridge_edges = matrix_patches_mean==10000;
        X_hard_ridge_edge = Xmesh_patch_regions(ix_hard_ridge_edges);
        Y_hard_ridge_edge = Ymesh_patch_regions(ix_hard_ridge_edges);



    
        for kk=1:length(Ygrid_vector)
            for kkk=1:n_windows

                P1=[Xgrid_vector(kk,:) Ygrid_vector(kk,:)];
                P2= window_xy(kkk,:);
                % define equation line
                ax = (P2(1,2)-P1(1,2))/(P2(1,1)-P1(1,1));
                bx = -1*ax*P1(1,1)+P1(1,2);
                ay =  (P2(1,1)-P1(1,1))/(P2(1,2)-P1(1,2));
                by = -1*ay*P1(1,2)+P1(1,1);
                x_sampling = [min([P1(1) P2(1)]) : max([P1(1) P2(1)])];
                y_sampling = [min([P1(2) P2(2)]) : max([P1(2) P2(2)])];
                y_line = (ax*x_sampling+bx);
                x_line = (ay*y_sampling+by);

                % limit the number of ridge intersections by a bounding box the two line points
                ix_temp =  (X_hard_ridge_edge >= min([P1(1) P2(1)])-2 & X_hard_ridge_edge <= max([P1(1) P2(1)])+2 & Y_hard_ridge_edge>= min([P1(2) P2(2)])-2 & Y_hard_ridge_edge<= max([P1(2) P2(2)])+2);

                 if sum(ix_temp)>0
                    % computing the residuals between the line and the ridge location
                    % x-direction sampling
                    y_line_res = reshape(repmat(y_line,length(Y_hard_ridge_edge(ix_temp)),1) - repmat(Y_hard_ridge_edge(ix_temp),1,length(y_line)),[],1);
                    x_sampling_res = reshape(repmat(x_sampling,length(X_hard_ridge_edge(ix_temp)),1) - repmat(X_hard_ridge_edge(ix_temp),1,length(x_sampling)),[],1);
                    XY_line_res = (abs([x_sampling_res y_line_res ]));
                    % x-direction sampling
                    x_line_res = reshape(repmat(x_line,length(X_hard_ridge_edge(ix_temp)),1) - repmat(X_hard_ridge_edge(ix_temp),1,length(x_line)),[],1);
                    y_sampling_res = reshape(repmat(y_sampling,length(Y_hard_ridge_edge(ix_temp)),1) - repmat(Y_hard_ridge_edge(ix_temp),1,length(y_sampling)),[],1);
                    YX_line_res = (abs([x_line_res y_sampling_res ]));

                    ix_ridge = sum(XY_line_res(:,1)<1 & XY_line_res(:,2)<1 );
                    ix_ridge = ix_ridge + sum(YX_line_res(:,1)<1 & YX_line_res(:,2)<1);
                else
                    ix_ridge = 0; 
                 end



                if ix_ridge>0
                    Distance_grid(kk,kkk)=NaN;
    %                 figure
    %                 plot(X_hard_ridge_edge,Y_hard_ridge_edge,'k.')
    %                 hold on
    %                 plot([P1(1) P2(1)],[P1(2) P2(2)],'r.-')
    %                 hold on
    %                 plot(X_hard_ridge_edge(ix_temp),Y_hard_ridge_edge(ix_temp),'go')
    % 
    %                 
    %                 axis equal
    %                 axis tight
    %                 keyboard
                end         
            end
        end
        fprintf('Computing the regional weights, Done ... \n')
        
        
        % find those poitns that were on the ridge edge and did not have a
        % window allocated. Set it to the closed window
        ix = find(sum(isnan(Distance_grid),2)==n_windows);
        if ~isempty(ix)
            % Computing weights based on distance from the points to the windows centers
            X_points = repmat(Xgrid_vector(ix,1),1,n_windows);
            Y_points = repmat(Ygrid_vector(ix,1),1,n_windows);
            X_windows = repmat(window_xy(:,1)',length(ix),1);
            Y_windows = repmat(window_xy(:,2)',length(ix),1);
            Distance_temp = sqrt((X_points-X_windows).^2+(Y_points-Y_windows).^2);
            clear X_points Y_points X_windows Y_windows

            Distance_temp_min = (min(Distance_temp'))';
            for kkk=1:length(Distance_temp_min)
                ix_closest_window = find(Distance_temp_min(kkk)==Distance_temp(kkk,:));
                Distance_grid(ix(kkk),ix_closest_window)=Distance_temp_min(kkk);
            end
        end


        if flag_figure_checks==1
            for kkk=1:n_windows
                figure('name',['Window ' num2str(kkk) ' distance to points (NaN is no influence due to ridge blockage)'])
                scatter3(Xgrid_vector(:,1),Ygrid_vector(:,1),Distance_grid(:,kkk),3,Distance_grid(:,kkk),'filled')
                view(0,90)
                axis equal
                axis tight
                title(['Window ' num2str(kkk) ' distance to points (NaN is no influence due to ridge blockage)'] )
                xlim([min(xy_local(:,1)) max(xy_local(:,1))]);
                ylim([min(xy_local(:,2)) max(xy_local(:,2))]);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%
        
        
        
        
        
    end
        
        
        
        
        
        
        
        
        
        
        
        
% % %         
% % %         
% % %         
% % %         % see if there are patches that are separated by windows.
% % %         ix_hard_ridge_edges = matrix_patches_mean==10000;
% % %         X_hard_ridge_edge = Xmesh_patch_regions(ix_hard_ridge_edges);
% % %         Y_hard_ridge_edge = Ymesh_patch_regions(ix_hard_ridge_edges);
% % % 
% % % 
% % %         fprintf('Computing the regional weights, this can take some time ... \n')
% % %         for kk=1:n_points
% % %             for kkk=1:n_windows
% % % 
% % %                 P1=xy_local(kk,:);
% % %                 P2= window_xy(kkk,:);
% % %                 % define equation line
% % %                 ax = (P2(1,2)-P1(1,2))/(P2(1,1)-P1(1,1));
% % %                 bx = -1*ax*P1(1,1)+P1(1,2);
% % %                 ay =  (P2(1,1)-P1(1,1))/(P2(1,2)-P1(1,2));
% % %                 by = -1*ay*P1(1,2)+P1(1,1);
% % %                 x_sampling = [min([P1(1) P2(1)]) : max([P1(1) P2(1)])];
% % %                 y_sampling = [min([P1(2) P2(2)]) : max([P1(2) P2(2)])];
% % %                 y_line = (ax*x_sampling+bx);
% % %                 x_line = (ay*y_sampling+by);
% % % 
% % %                 % limit the number of ridge intersections by a bounding box the two line points
% % %                 ix_temp =  (X_hard_ridge_edge >= min([P1(1) P2(1)])-2 & X_hard_ridge_edge <= max([P1(1) P2(1)])+2 & Y_hard_ridge_edge>= min([P1(2) P2(2)])-2 & Y_hard_ridge_edge<= max([P1(2) P2(2)])+2);
% % % 
% % %                  if sum(ix_temp)>0
% % %                     % computing the residuals between the line and the ridge location
% % %                     % x-direction sampling
% % %                     y_line_res = reshape(repmat(y_line,length(Y_hard_ridge_edge(ix_temp)),1) - repmat(Y_hard_ridge_edge(ix_temp),1,length(y_line)),[],1);
% % %                     x_sampling_res = reshape(repmat(x_sampling,length(X_hard_ridge_edge(ix_temp)),1) - repmat(X_hard_ridge_edge(ix_temp),1,length(x_sampling)),[],1);
% % %                     XY_line_res = (abs([x_sampling_res y_line_res ]));
% % %                     % x-direction sampling
% % %                     x_line_res = reshape(repmat(x_line,length(X_hard_ridge_edge(ix_temp)),1) - repmat(X_hard_ridge_edge(ix_temp),1,length(x_line)),[],1);
% % %                     y_sampling_res = reshape(repmat(y_sampling,length(Y_hard_ridge_edge(ix_temp)),1) - repmat(Y_hard_ridge_edge(ix_temp),1,length(y_sampling)),[],1);
% % %                     YX_line_res = (abs([x_line_res y_sampling_res ]));
% % % 
% % %                     ix_ridge = sum(XY_line_res(:,1)<1 & XY_line_res(:,2)<1 );
% % %                     ix_ridge = ix_ridge + sum(YX_line_res(:,1)<1 & YX_line_res(:,2)<1);
% % %                 else
% % %                     ix_ridge = 0; 
% % %                  end
% % % 
% % % 
% % % 
% % %                 if ix_ridge>0
% % %                     Distance(kk,kkk)=NaN;
% % %     %                 figure
% % %     %                 plot(X_hard_ridge_edge,Y_hard_ridge_edge,'k.')
% % %     %                 hold on
% % %     %                 plot([P1(1) P2(1)],[P1(2) P2(2)],'r.-')
% % %     %                 hold on
% % %     %                 plot(X_hard_ridge_edge(ix_temp),Y_hard_ridge_edge(ix_temp),'go')
% % %     % 
% % %     %                 
% % %     %                 axis equal
% % %     %                 axis tight
% % %     %                 keyboard
% % %                 end         
% % %             end
% % %         end
% % %         fprintf('Computing the regional weights, Done ... \n')
% % % 
% % %     end     
% % %     
% % %     % find those poitns that were on the ridge edge and did not have a
% % %     % window allocated. Set it to the closed window
% % %     ix = find(sum(isnan(Distance),2)==n_windows);
% % %     if ~isempty(ix)
% % %         % Computing weights based on distance from the points to the windows centers
% % %         X_points = repmat( xy_local(ix,1),1,n_windows);
% % %         Y_points = repmat(xy_local(ix,2),1,n_windows);
% % %         X_windows = repmat(window_xy(:,1)',length(ix),1);
% % %         Y_windows = repmat(window_xy(:,2)',length(ix),1);
% % %         Distance_temp = sqrt((X_points-X_windows).^2+(Y_points-Y_windows).^2);
% % %         clear X_points Y_points X_windows Y_windows
% % % 
% % %         Distance_temp_min = (min(Distance_temp'))';
% % %         for kkk=1:length(Distance_temp_min)
% % %             ix_closest_window = find(Distance_temp_min(kkk)==Distance_temp(kkk,:));
% % %             Distance(ix(kkk),ix_closest_window)=Distance_temp_min(kkk);
% % %         end
% % %     end
% % %     
% % % 
% % %     if flag_figure_checks==1
% % %         for kkk=1:n_windows
% % %             figure('name',['Window ' num2str(kkk) ' distance to points (NaN is no influence due to ridge blockage)'])
% % %             scatter3(xy_local(:,1),xy_local (:,2),Distance(:,kkk),3,Distance(:,kkk),'filled')
% % %             view(0,90)
% % %             axis equal
% % %             axis tight
% % %             title(['Window ' num2str(kkk) ' distance to points (NaN is no influence due to ridge blockage)'] )
% % %             xlim([min(xy_local(:,1)) max(xy_local(:,1))]);
% % %             ylim([min(xy_local(:,2)) max(xy_local(:,2))]);
% % %         end
% % %     end
    
   
    % Compute the gaussian weight based on distance
    w_d = normpdf(Distance_grid,0,patch_c2c);                               % [n_points x n_windows]
    
    % Put the weight of the patches excluded due to mountains to zero.
    w_d(isnan(w_d))=0;

    % storing other variables specific for the other computation method
    window_xy_extra = [];
    ix_edges_window =[];
    ridges_flag = 'y';
    XYgrid_vector = [Xgrid_vector Ygrid_vector];
    n_grid_average = ceil(2*patch_c2c./grid_resolution);

    % getting the variables in a structure to load on the next run
    iterate.window_ix = window_ix;
    iterate.window_xy = window_xy;
    iterate.window_xy_extra = window_xy_extra;
    iterate.ix_edges_window = ix_edges_window;
    iterate.w_d=w_d;
    iterate.ridges_flag=ridges_flag;
    iterate.patch_regions=patch_regions;
    iterate.Ymesh_patch_regions = Ymesh_patch_regions;
    iterate.Xmesh_patch_regions = Xmesh_patch_regions;
    iterate.XYgrid_vector = XYgrid_vector;
    iterate.grid_resolution = grid_resolution;
    iterate.n_grid_average = n_grid_average;
    iterate.dimensions_grid = dimensions_grid;
    clear Xgrid_vector Ygrid_vector

else
    % Loading the existing data
    % When using the ridges approach some of the fields are empty.
    % All patches wil have an equal weight in the ridges approach
    % All points within a patch will have the patche value. 
    % But a mean filter will be used to smooth out the data
    window_ix = iterate.window_ix;
    window_xy = iterate.window_xy;
    window_xy_extra = iterate.window_xy_extra;
    ix_edges_window = iterate.ix_edges_window;
    w_d = iterate.w_d; 
    ridges_flag = iterate.ridges_flag;
    patch_regions = iterate.patch_regions;
    Ymesh_patch_regions = iterate.Ymesh_patch_regions;
    Xmesh_patch_regions = iterate.Xmesh_patch_regions ;
    
    if  strcmpi(powerlaw_ridge_constraint,'y')
        XYgrid_vector = iterate.XYgrid_vector;
        grid_resolution = iterate.grid_resolution;
        n_grid_average = iterate.n_grid_average;
        dimensions_grid = iterate.dimensions_grid;
    end

    
end
clear kk k


%% Estimating the linear relation for each window.
% do it for each frequency band and keep only those which are consistent
n_windows = length(window_ix);
n_datasets = size(data_local1,2);


% initialisation of the variables
slope_window = NaN([n_windows 1]);
% looping over all the windows
for k=1:n_windows
    if n_windows==1 && strcmp(plane_mode,'y')  % estimate a planar K using all points in a single window
        ix = window_ix{k};
        % number of points in the window
        n_points_window = length(ix);
        
  
        
        for kk=1:n_datasets        
            
            if sum(isnan(data_local1(ix,kk)))>1
                fprintf('Modify the inversion such no nans are pressent, this is for crop out region! \n')
                keyboard 
            end

            
            scaling = 1./mean(abs(data_local1(ix,kk)));
            % set up the design matrix 
            A = [data_local1(ix,kk).*scaling.*xy_local(:,1) data_local1(ix,kk).*scaling.*xy_local(:,2) data_local1(ix,kk).*scaling ones([n_points_window 1])];
            % doing the (equal weight) least square inversion
            coeff = lscov(A,data_local2(ix,kk));
            % correct for the introduced scaling
            slope_window_band(kk,1) = coeff(1).*scaling;
            slope_window_band(kk,2) = coeff(2).*scaling;
            slope_window_band(kk,3) = coeff(3).*scaling;
            offset_window_band(kk) = coeff(4);
            
            % popagation of the variances for the estimated coefficients
            Q_coeff = inv(A'*A);
            % correct for the introduced scaling
            slope_window_band_std(kk,1) = sqrt(Q_coeff(1,1)).*scaling;
            slope_window_band_std(kk,2) = sqrt(Q_coeff(2,2)).*scaling;
            slope_window_band_std(kk,3) = sqrt(Q_coeff(3,3)).*scaling;
            offset_window_band_std(kk) = sqrt(Q_coeff(2,2));      
        end
        clear kk

        % Checking if the troposheric delay needs to be computed for each bandfilter
        if strcmp(powerlaw_all_bands,'y')
            % keep the information of all the band stored for latter on
            slope_window_band_matrix(:,k) = slope_window_band;              % matrix [n_bands 3 components of plane]
            slope_window_band_std_matrix(:,k) = slope_window_band_std;      % matrix [n_bands 3 components of plane]
        end
        
        % Estimating the window slope value by the weighted mean of the consistent bands.
        % deformation contaminated bands are rejected by outlier testing.
        % outlier rejection based on the bands for each of the plane coefficients
        [slope_mean_x,slope_window_std_mean_x,obs_kept] = w_test_mean(slope_window_band(:,1),slope_window_band_std(:,1),[]);
        slope_window(k,1) = slope_mean_x;
        slope_window_std(k,1) = slope_window_std_mean_x;
        [slope_mean_y,slope_window_std_mean_y,obs_kept] = w_test_mean(slope_window_band(:,2),slope_window_band_std(:,2),[]);
        slope_window(k,2) = slope_mean_y;
        slope_window_std(k,2) = slope_window_std_mean_y;
        [slope_mean,slope_window_std_mean,obs_kept] = w_test_mean(slope_window_band(:,3),slope_window_band_std(:,3),[]);
        slope_window(k,3) = slope_mean;
        slope_window_std(k,3) = slope_window_std_mean;
        clear slope_mean_x slope_window_std_mean_x obs_kept slope_mean_y slope_window_std_mean_y slope_mean slope_window_std_mean
        

        
    else % estimate a single value for K in each subwindow
        ix_original = window_ix{k};
        % number of points in the window
        n_points_window_original = length(ix_original);

        
        
        % looping over the different bandfiltered datasets
        % initialisation of the variables
        slope_window_band = NaN([n_datasets 1]);
        slope_window_band_std = NaN([n_datasets 1]);
        slope_window_band_std_bootstrap = NaN([n_datasets 1]);
        slope_window_band_std_error_prop = NaN([n_datasets 1]);
        offset_window_band = NaN([n_datasets 1]);
        offset_window_band_std = NaN([n_datasets 1]);
        
        
        for kk=1:n_datasets
            
            % removing the nan values for each data set, which could be
            % different for each interferogram, or could be a mutially
            % cropped region
            ix_bad = ((isnan(data_local1(ix_original,kk)) + isnan(data_local2(ix_original,kk)))>=1);
%             ix = ix_original(~isnan(data_local1(ix_original,kk)));
            ix = ix_original(~ix_bad);
            n_points_window = length(ix);
            
            if n_points_window>=10
                % use mean to renormalize such the matrix does not become rank
                % deficient for extreme values
                scaling = 1./mean(abs(data_local1(ix,kk)));
                
                 
                % set up the design matrix 
                A = [data_local1(ix,kk).*scaling ones([n_points_window 1])];
                % doing the (equal weight) least square inversion
                coeff = lscov(A,data_local2(ix,kk));
                
                if sum(isnan(coeff))>0
                    keyboard
                end
                % correct for the introduced scaling
                slope_window_band(kk) = coeff(1).*scaling;
                offset_window_band(kk) = coeff(2);
                clear coeff 


             

                % bootstrapping for error estimation    
                [coeff_std,coeff_vector]= aps_powerlaw_bootstrap(A,data_local2(ix,kk),n_boot_runs);
                slope_window_band_std_bootstrap(kk)=coeff_std(1).*scaling;

         
                % popagation of the variances for the estimated coefficients
                Q_coeff = inv(A'*A);
                % correct for the introduced scaling
                slope_window_band_std_error_prop(kk) = sqrt(Q_coeff(1,1)).*scaling;          
                clear Q_coeff scaling A coeff
                
            else
                % not enough points to estimate the relationship
                slope_window_band_std_bootstrap(kk) = NaN;
                slope_window_band_std_error_prop(kk) = NaN;

            end
        end
        clear kk ix
        
        
        % ifg you want to change the std you should do it here:
        if bootstrap_flag==1
           slope_window_band_std = slope_window_band_std_bootstrap;      % matrix [n_bands n_windows]  
        else
           slope_window_band_std = slope_window_band_std_error_prop;      % matrix [n_bands n_windows]              
        end
           
        % Checking if the troposheric delay needs to be computed for each bandfilter
        if strcmp(powerlaw_all_bands,'y')
            % keep the information of all the band stored for latter on
            slope_window_band_matrix(:,k) = slope_window_band;              % matrix [n_bands n_windows]
            slope_window_band_std_matrix(:,k) = slope_window_band_std;      % matrix [n_bands n_windows]
            
        end



        % Estimating the window slope value by the weighted mean of the consistent bands.
        % deformation contaminated bands are rejected by outlier testing.
        [slope_mean,slope_window_std_mean,obs_kept] = w_test_mean(slope_window_band,slope_window_band_std,[]);
        slope_window(k,1) = slope_mean;
        slope_window_std(k,1) = slope_window_std_mean;
        

        if save_flag_outlier_windows==1  && strcmpi(powerlaw_ridge_constraint,'n')
           temp = [1:n_datasets]';
           edge_flag = find(k==ix_edges_window);
           if isempty(edge_flag)~=1
               str_extra = ' (edge window)';
           else
              str_extra = ''; 
           end
           hfig = figure('name',['outlier rejection for window' num2str(k) str_extra]);
           errorbar(slope_window_band,slope_window_band_std,'r.')
           hold on
           plot(temp(obs_kept),slope_window_band(obs_kept),'go')
           hold on
           plot([1 n_datasets],[slope_mean slope_mean],'k-')
           legend('obs','obs kept','mean estimate')

           if exist(save_path_outlier_window,'dir')~=7
                mkdir(save_path_outlier_window)
           end
           title({['window ' num2str(k)],' '})
           xlabel(['Band 1 <----     -----> Band n'])


           if isempty(dimension_filter)~=1
                Ax1 = gca;
                ylims = get(Ax1,'ylim');

                ixplot = find(dimension_filter==1);
                if isempty(ixplot)~=1
                    filter_transition = ixplot(1);
                    hold on
                    plot([filter_transition filter_transition],ylims,'k--')
                    legend('obs','obs kept','mean estimate','2D - 1D transition')

                end
                clear ixplot
           end
           set(hfig,'PaperPositionMode','auto')
           print(hfig,'-dpng','-r150',[save_path_outlier_window filesep 'outlier_rej_window_' num2str(k)])
           close(hfig)
           clear temp hfig;
        end
        clear slope_mean obs_kept
    end
    fprintf(['Completed window ' num2str(k) ' out of ' num2str(n_windows) '\n'])
end


perc_cons_K=NaN;
if n_windows==1 && strcmp(plane_mode,'y')  
    % Checking if the troposheric delay needs to be computed for each bandfilter
    if strcmp(powerlaw_all_bands,'y')
        for kk=1:n_datasets
            
            if sum(isnan(slope_window_band_matrix(:,kk)))>1
                fprintf('Modify the inversion such no nans are pressent, this is for crop out region! \n')
                keyboard 
            end
        
            
            
    
           % Compute slope for each PS point and the long wavelength tropospheric signal for each ifgs
            slope_band = [xy_local ones([size(xy_local(:,1),1) 1])]*slope_window_band_matrix(:,kk);     % note that slope_window_band_matrix is defined slightly different for 1 or n mindows
           
            % compute delay from it
            ph_tropo_powerlaw_band(:,kk) = slope_band.*(h0(ifg_number)*1000-hgt).^(alpha(ifg_number));
            
            slope_local_band(:,kk) = slope_band;
            clear slope_band
        end
    else
        ph_tropo_powerlaw_band = [];
        slope_local_band = [];
    end
    % computing K at each PS location using the planar estimate
    slope_local = [xy_local ones([size(xy_local(:,1),1) 1])]*slope_window';
    
    fprintf('Plane mode \n')
    
else
    % Checking if the troposheric delay needs to be computed for each bandfilter
    if strcmp(powerlaw_all_bands,'y')
        for kk=1:n_datasets
                
%             keyboard
            
                if sum(isnan(slope_window_band_std_matrix(kk,:)))>1
                    fprintf('Modify the inversion such no nans are pressent, this is for crop out region! \n')
                    keyboard 
                end
            
            
            
            
                w_std = 1./slope_window_band_std_matrix(kk,:);                  % Large std -> low weight
                w_std = repmat(w_std,size(w_d,1),1);
                w_std((w_d==0))=0;
                w_std = w_std./repmat(sum(w_std,2),1,n_windows);

                % combine the weights on distance and standard deviation of the windows together:
                w = w_d.*w_std;%repmat(w_std,n_points,1);
                w = w./repmat(sum(w,2),1,n_windows);
                clear w_std

                if flag_figure_checks==1  && strcmpi(powerlaw_ridge_constraint,'n')
                    figure
                    scatter3(iterate.window_xy(:,1),iterate.window_xy(:,2),w(1000,:),25,w(1000,:),'filled')
                    view(0,90)
                    axis equal
                    axis tight
                    colorbar
                end

                % For the ridges approach the data is kept on a regular
                % grid as at the end a smoothing is applied, after which
                % the interpolation to the ps lcoations is done.
                % For the windows approach the weights are at onces based
                % on those for the PS point locations
                if strcmpi(powerlaw_ridge_constraint,'y')
                    
                    slope_band_grid = w*slope_window_band_matrix(kk,:)';
                    % putting it back in a regular grid.
                    temp_matrix = reshape(slope_band_grid,dimensions_grid(1),[]);
                    
                    % mean filter accros hard ridges
                    H = fspecial('average',n_grid_average );
                    temp_matrix_mean = imfilter(temp_matrix,H,'replicate');
                    
                    slope_band_grid_mean =reshape(temp_matrix_mean,[],1);

                    if flag_figure_checks==1
                        figure
                        scatter3(XYgrid_vector(:,1),XYgrid_vector(:,2),slope_band_grid,25,slope_band_grid,'filled')
                        view(0,90)
                        axis equal
                        axis tight
                        colorbar
                        
                        figure
                        scatter3(XYgrid_vector(:,1),XYgrid_vector(:,2),slope_band_grid_mean,25,slope_band_grid_mean,'filled')
                        view(0,90)
                        axis equal
                        axis tight
                        colorbar
                        
                    end
                    % perform smoothing operation
                    slope_band = griddata(XYgrid_vector(:,1),XYgrid_vector(:,2),slope_band_grid_mean,xy_local(:,1),xy_local(:,2),'linear');

                else
                    % Compute slope, c for each PS point and the long wavelength tropospheric signal for each ifgs
                    slope_band = w*slope_window_band_matrix(kk,:)';       % excluding the master
                end
                clear  w
           
                % compute delay from it
                ph_tropo_powerlaw_band(:,kk) = slope_band.*(h0(ifg_number)*1000-hgt).^(alpha(ifg_number));
                slope_local_band(:,kk)=slope_band;
                clear slope_band
        end
        clear kk slope_window_band_std_matrix  slope_window_band_matrix
    else
        ph_tropo_powerlaw_band = [];
        slope_local_band = [];
    end
    
%     keyboard

    %% Computing the delay based on the slope estimation from the good bands     
    % The weights based on the std of the windows
    w_std = 1./slope_window_std';                  % Large std -> low weight
    w_std =  repmat(w_std,size(w_d,1),1);
    w_std((w_d==0))=0;
    w_std = w_std./repmat(sum(w_std,2),1,n_windows);

    % combine the weights on distance and standard deviation of the windows together:
    w = w_d.*w_std;
    w = w./repmat(sum(w,2),1,n_windows);
    clear w_std

    % Compute slope, c for each PS point and the long wavelength tropospheric signal for each ifgs
    if strcmpi(powerlaw_ridge_constraint,'y')
        slope_grid = w*slope_window;
        
        ix_pos = slope_window>=0;
        perc_cons_K = floor((abs(sum(ix_pos)./length(slope_window)-0.5)+0.5)*100);
        fprintf([num2str(perc_cons_K) ' perc of K are consistent in space \n ']);    
        
        % putting it back in a regular grid.
        temp_matrix = reshape(slope_grid,dimensions_grid(1),[]);

        % mean filter accros hard ridges
        H = fspecial('average', n_grid_average);
        temp_matrix_mean = imfilter(temp_matrix,H,'replicate');

        slope_grid_mean =reshape(temp_matrix_mean,[],1);

        % perform smoothing operation
        slope_local = griddata(XYgrid_vector(:,1),XYgrid_vector(:,2),slope_grid_mean,xy_local(:,1),xy_local(:,2),'linear');
    else
        perc_cons_K=nan;
        slope_local = w*slope_window;
    end
    clear w w_d
    
    
    if flag_figure_checks==1 && debug_fig==1
       figure
       scatter3(xy_local(:,1),xy_local(:,2),slope_local,3,slope_local,'filled')
       hold on
       scatter3(window_xy(:,1),window_xy(:,2),slope_window,50,slope_window,'filled')
       view(0,90)
       axis equal
       axis tight
       xlimits = get(gca,'xlim');
       ylimits = get(gca,'ylim');     
       title('slope')

       figure
       scatter3(window_xy(:,1),window_xy(:,2),slope_window_std,50,slope_window_std,'filled')
       view(0,90)
       axis equal
       axis tight
       title('std')
       set(gca,'xlim',xlimits)
       set(gca,'ylim',ylimits)

    end
    
    fprintf('Regular mode \n')
end

% compute delay from it
ph_tropo_powerlaw = slope_local.*(h0(ifg_number)*1000-hgt).^(alpha(ifg_number));
if ~isempty(ix_phase_nan_interferogram)
    ph_tropo_powerlaw(ix_phase_nan_interferogram,:)=NaN;
end
