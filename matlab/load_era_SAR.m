function [xyz_input,xyz_output] = load_era_SAR(filename,xy_out_grid,overwrite_flag) 
% [xyz] = load_era_SAR(filename,crop_filename,geocoord_flag) 
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
%                       assumed. if xyz_output is inf, then the data has
%                       been stored and not given as output. This is
%                       default option for grd files, as they can become
%                       large.

% By David Bekaert - University of Leeds
%
% Modifications:
% DB    10/06/2015      Fix in case the lonlat observations are a rectangle
% DB    24/11/2015      Include chunking of the data for large datasets
% DB    24/11/2015      Add support for grd files
% DB    30/11/2015      Remove the polygon cropping option, introduce
%                       overwrite flag overwrite_flag instead. Include the
%                       option to use MEX file to write out the data. This
%                       is much faster than standard matlab, but requires
%                       correct compiling of the mex code.


plot_flag = 0;

% options of handlign larger data.
% 1) use of .grd files and gmt. In future the .grd approach will become the
% default processing in aps_era_SAR and aps_wrf_SAR.
% this however requires the data grid to be written once to a txt file,
% which with regular matlab can take long time. Therefore there is the
% option to use a mex variant that is much faster, but needs to be compiled
% first. By default this will be tried, and reversed to conventional
% data writing if it fails.

% 2) use of .xyz datafiles and matlab. This is backward compatible with old
% TRAIN processed data. For large datasets this will becomes challenging,
% instead you can chunk you data to make it less memory intensive. By
% default this option is turned on. If set to 'n' matlab will interpolate
% all data at once.
chunk_flag = 'y';       % when yes use the chunk approach to interpolate to the ps locations
chunk_bin_size = 0.5;   % size in deg of the lonlat bins used for chunkign the data

if nargin <1
    error('myApp:argChk', ['Not enough input arguments...\n'])    
end
if nargin <2 || isempty(xy_out_grid)
    interpolate_flag = 0;    
else
    interpolate_flag = 1;    
end
if nargin <3 || isempty(overwrite_flag)
    overwrite_flag = 1;    
end

    


%% checking the file type if its grd or not
file_type_grd = 'n';
[file_path,temp,file_ext] = fileparts(filename);
% check if this is a grd file already
if strcmpi(file_ext,'.grd')
    fprintf(['Specified aps correction is a .grd file\n'])
    file_type_grd = 'y';
    xyz_input = inf;
    
    [gmt5_above, gmt_version]=get_gmt_version();
    error('This is currently under development')    
    
else
    %% Loading of processed data files
    fid = fopen(filename,'r');
    data_vector = fread(fid,'double');
    fclose(fid);
    
    % reshaping into the right 3 column matrix
    xyz_input = reshape(data_vector,3,[])';
    clear data_vector
end



%% Crop the dataset when usign a scatter approach in matlab
% When doing the interpolation, crop the data first to the maximum required extend.
% overwrite any specific crop file then.
if ~strcmpi(file_ext,'.grd')
    if interpolate_flag ==1
       crop_flag =1;
       xy_min = min(xy_out_grid);
       xy_max = max(xy_out_grid);
       xy_min = xy_min-0.01;
       xy_max = xy_max+0.01;

       poly.xy = [xy_min(1)   xy_min(2);
                  xy_min(1)   xy_max(2);
                  xy_max(1)   xy_max(2);
                  xy_max(1)   xy_min(2)];
        % adapt the grid when need to have the same coordinates as the crop
        if xy_min(1)<0 & min(xyz_input(:,1))>0
            xyz_input(:,1)=xyz_input(:,1)-360;
        end


        % do the actual cropping
        ix = inpolygon( xyz_input(:,1) , xyz_input(:,2) , poly.xy(:,1) ,poly.xy(:,2) ) ;
        % including those points on the polygon edge as well
        ix(ix>1)=1;

        % keeping only those data within the cropaps
        xyz_output = xyz_input(ix,:);

        clear ix
    end
end


%% When requested interpolate the data
if interpolate_flag==1
    if strcmpi(file_ext,'.grd')
 
        
        % writing out the lonlat coordiantes in an ASCI table for GMT to use latter on
        lonlat_file = [file_path filesep '..' filesep 'lonlat_temp.txt'];
           
        
        % check if the lonlat file needs to be made
        generate_file = 'n';
        if  overwrite_flag==0 && exist(lonlat_file,'file')~=2
            generate_file='y';
        end
        if overwrite_flag==1
            generate_file='y';
        end
                
        % This file writign is intensive for high resolution data files.
        % The MEX file uses c-code and is much faster, but needs to be
        % compiled for your system. By default the slow matlab approach is
        % used. But once compiled you can change this flag to use the
        % faster code.
        if strcmpi(generate_file,'y')
            try 
                dumptofile(xy_out_grid, lonlat_file,' '); 
            catch 
                dlmwrite(lonlat_file,xy_out_grid);
            end
        end
        
        % generating the ifg filename
        [SAR_path,SAR_filename,temp] = fileparts(filename);
        save_name = [SAR_path filesep SAR_filename '_SARll'];
        % calling GMT for the interpolation to the ifg lonlat grid
        if strcmpi(gmt5_above,'y')
           commandstr = ['grdtrack ' lonlat_file ' -G' filename '  -fg -N  -Z > ' save_name]
        else
           commandstr = ['grdtrack ' lonlat_file ' -G' filename '  -fg -Qn  -Z > ' save_name]
            
        end
        aps_systemcall(commandstr);

        % turn the plot grid off as this is not supported for this case (large data amount)
        plot_flag = 0;

        % give the output grid a 'inf' flag such it can be used to trigger the code
        % that this is data saved.
        xyz_output = [save_name];
        
   
    elseif strcmpi(chunk_flag,'y')
        
       % interpolate in chunks when requested
       xyz_output_temp = xyz_output;
       clear xyz_output
       xyz_output = NaN([size(xy_out_grid,1) 3]);

       % define a grid which will be used to loop over 
       xy_min = floor(min(xy_out_grid));
       xy_max = ceil(max(xy_out_grid));
       
       % the grid step size depends on the chunck size
       lon_grid = [xy_min(1):chunk_bin_size:xy_max(1)-chunk_bin_size];
       lat_grid = [xy_min(2):chunk_bin_size:xy_max(2)-chunk_bin_size];
       [lon_grid,lat_grid] = meshgrid(lon_grid,lat_grid);
       lon_grid = reshape(lon_grid,[],1);
       lat_grid = reshape(lat_grid,[],1);

       % start the loop
       tic
       fprintf(['Buffering in ' num2str(chunk_bin_size) ' deg chunks... \n'])
       for k=1:length(lat_grid)
           ix_temp = (xy_out_grid(:,1)>=lon_grid(k) & xy_out_grid(:,1)<=lon_grid(k)+chunk_bin_size & xy_out_grid(:,2)>=lat_grid(k) & xy_out_grid(:,2)<=lat_grid(k)+chunk_bin_size);
           if sum(ix_temp)>0
               ix_temp_full = (xyz_output_temp(:,1)>=lon_grid(k)-chunk_bin_size/2 & xyz_output_temp(:,1)<=lon_grid(k)+chunk_bin_size+chunk_bin_size/2 & xyz_output_temp(:,2)>=lat_grid(k)-chunk_bin_size/2 & xyz_output_temp(:,2)<=lat_grid(k)+chunk_bin_size+chunk_bin_size/2);
               z_output = griddata(xyz_output_temp(ix_temp_full,1),xyz_output_temp(ix_temp_full,2),xyz_output_temp(ix_temp_full,3),xy_out_grid(ix_temp,1),xy_out_grid(ix_temp,2),'linear');
               xyz_output(ix_temp,:) = [xy_out_grid(ix_temp,:) z_output] ;
           end
           fprintf([num2str(k) '/' num2str(length(lat_grid)) '\n']);
       end
       toc
        
        
    else
        xyz_output_temp = xyz_output;
        clear xyz_output
        z_output = griddata(xyz_output_temp(:,1),xyz_output_temp(:,2),xyz_output_temp(:,3),xy_out_grid(:,1),xy_out_grid(:,2),'linear');
        xyz_output = [xy_out_grid z_output] ;
    end
end


%% Empty output variable in case no of the optional commands are selected
if interpolate_flag==0
    xyz_output = [];
else
    if plot_flag==1
        % plotting the cropped and or interpolated result
        figure('name','Cropped and or interpolated dataset')
        scatter3(xyz_output(:,1),xyz_output(:,2),xyz_output(:,3),3,xyz_output(:,3),'filled')
        view(0,90)
        axis equal
        axis tight
        colorbar
        box on
    end
end



if plot_flag==1
    % plotting the cropped and or interpolated result
    figure('name','Input dataset')
    scatter3(xyz_input(:,1),xyz_input(:,2),xyz_input(:,3),3,xyz_input(:,3),'filled')
    view(0,90)
    axis equal
    axis tight
    colorbar
    box on
end


