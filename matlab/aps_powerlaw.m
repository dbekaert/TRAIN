function aps_powerlaw(start_step,end_step,save_path)
% aps_powerlaw(start_step,end_step)
% Program to compute the powerlaw tropsoheric delays. 
%   STEP 1 = Estimate powerlaw delay and maximum height from sounding data
%            (when this is available), or use defaults or user set values
%   STEP 2 = Rotate dataset to minimize interpolation effects, Powerlaw scaling and interpolate to regular grid 
%   STEP 3 = 1D/2D bandfiltering in the frequency domain and interpolate 
%            back to local grid 
%   STEP 4 = Spatial local estimation of the powerlaw scaling coefficient,
%            estimate final values from reliable frequency band, extrapolate 
%            to all point locations and compute powerlaw tropospheric delay. 
% 
%   To change the processing parameters use setparm_aps and getparm_aps
%
% **** when using this estimation method please cite:
%       D.P.S. Bekaert, A.J. Hooper and T.J. Wright, A spatially-variable power-law tropospheric correction
% 	  	 technique for InSAR data, JGR, doi:10.1029/2014JB011558 *****
%
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
%
% modifications:
% 04/2013   DB:     Incorportate a parm_aps list with all the processing
%                   options and variables.
% 05/2013   DB:     Include powerlaw coefficient sensitivity analysis
% 05/2013   DB:     Allow for an interferogram based estimation
% 12/2013   DB:     Put the flag of ifg correction on prior to powerlaw computation
% 02/2014   DB:     Integrate the plotting with the running of the sounding data
% 02/2014   DB:     Include a mountain ridge code, user interaction needed
% 07/2014   DB:     Have the option to plot the ridges
% 07/2014   DB:     Include a cropping option
% 11/2014   DB:     Fix bug and include paper reference
% 03/2015   DB:     Ask for the reference technique in step 5
% 01/2016   DB:     Fix ll.lonlat, ph_uw, hgt for non-stamps data, replace
%                   save calls with aps_save, in case no -append flag is included
% 02/2016   DB:     Expand option 5 to plot with the unwrapped phase too.


if nargin<3
   save_path = [pwd filesep]; 
end

currdir = pwd;		% current working directory
plot_flag = 0;
hydro =1;
wet=1;

fprintf(['\n********************************************************************************************************** \n','D.P.S. Bekaert, A.J. Hooper and T.J. Wright,\n', 'A spatially-variable power-law tropospheric correction technique for InSAR data, JGR,  doi:10.1029/2014JB011558\n', '**********************************************************************************************************\n\n'])


%% PART 1: Define the powerlaw coefficents, from the user or by using sounding data.
sounding_data = getparm_aps('sounding_data',1);
stamps_processed = getparm_aps('stamps_processed',1);
powerlaw_ridge_constraint = getparm_aps('powerlaw_ridge_constraint',1);

if strcmp(stamps_processed,'y')
    load psver
else
    psver = 2; 
end

% file names of the output data
apsname = [save_path filesep 'tca' num2str(psver) '.mat'];
apssbname = [save_path filesep 'tca_sb' num2str(psver) '.mat'];
apsbandsname = [save_path filesep 'tca_bands' num2str(psver) '.mat'];
apsbandssbname = [save_path filesep 'tca_bands_sb' num2str(psver) '.mat'];


% saving part of the data in a subfolder
save_path = [save_path filesep 'aps_p'];
if exist(save_path,'dir')~=7
    mkdir(save_path);
end

if start_step==0 && strcmp(powerlaw_ridge_constraint,'y')
    
    mountain_ridge=[];
    if exist('tca_support.mat')==2
        load('tca_support.mat','powerlaw_ridges')
        if exist('powerlaw_ridges','var')==1
            mountain_ridge = powerlaw_ridges.mountain_ridge;
        end
    end
    if isempty(mountain_ridge)==1
        str_repr='y';
        str_vis = 'n';
    end
    
    % data has been computed before - reporcess? - visualize?
    if length(mountain_ridge)>0
        str_repr = '';
        while strcmpi(str_repr,'y')~=1 && strcmpi(str_repr,'n')~=1
            str_repr = input(['Ridges have been defined before, reprocess the data? [y/n] \n'],'s');
        end
        
        % no reprocessing, do you want to visualize it?
        if strcmpi(str_repr,'n')
            str_vis = '';
            while strcmpi(str_vis,'y')~=1 && strcmpi(str_vis,'n')~=1
                str_vis = input(['Do you want to visualize the data? [y/n] \n'],'s');
            end
        else
            str_vis = 'n';
        end
        
    end
    
    if strcmpi(str_repr,'y')
        [mountain_ridge] = aps_powerlaw_watershed;
    end
    if strcmpi(str_vis,'y')
        aps_support_plot(1);
    end
end

if start_step==1 && strcmp(sounding_data,'y')
    fprintf('\n\nStep 1: Powerlaw coefficients from the sounding data\n')
   
   
    % getting the variables from the parm_aps file
    look_angle = getparm_aps('look_angle',1);
    lambda = getparm_aps('lambda',1);
    sounding_h0= getparm_aps('sounding_h0',1);
    sounding_h_alpha_thres= getparm_aps('sounding_h_alpha_thres',1);
    sounding_start_date= getparm_aps('sounding_start_date',1);
    sounding_end_date= getparm_aps('sounding_end_date',1);
    sounding_time_stamp= getparm_aps('sounding_time_stamp',1);
    sounding_dir= getparm_aps('sounding_dir',1);
    sounding_error_promp= getparm_aps('sounding_error_promp',1);
    sounding_sensitivity= getparm_aps('sounding_sensitivity',1);
    n_months = getparm_aps('sounding_months',1);
    time_stamp = getparm_aps('sounding_time_stamp',1);
    sounding_ifg_dates = getparm_aps('sounding_ifg_dates',1);
    
    %%% check if the file is already exisiting - ask if reprocessing is
    %%% needed, if not load exisitng file and plot results 
    time_stamp_str = [];
    for k=1:size(time_stamp,1)
        if k>1
            time_stamp_str = [time_stamp_str '_' time_stamp(k,:)];
        else
            time_stamp_str = [time_stamp(k,:)];
        end
    end
    
    
    if strcmp(sounding_sensitivity,'y')
        fprintf('Sensitivity analyses of sounding data \n\n')
        
        % check first if the output file is already exisiting otherwize
        % ask if it needs to be reprocessed or visualised.
        
        if strcmp(sounding_ifg_dates,'y')
            fprintf('Coefficient estimation for each interferogram \n\n')
            if hydro==1 && wet==0
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_hydro_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
            elseif hydro==0 && wet==1
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_wet_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
            else
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
            end
            n_months =1;
        else
            fprintf('Coefficient estimation based on average \n\n')
            % putting the variables in the right set-up
            start_year = str2num(sounding_start_date(1:4));
            end_year = str2num(sounding_end_date(1:4));
            start_str = sounding_start_date(5:6);
            end_str = sounding_end_date(5:6);

            % the file name to be loaded
            if hydro==1 && wet==0
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_hydro_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
            elseif hydro==0 && wet==1
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_wet_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
            else
                save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
            end
        end


        if exist(save_name,'file')==2
            str = '';
            while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
                str = input(['This file already exist, do you want to re-process the data? [y/n] \n'],'s');
            end
            
            if strcmpi(str,'y')
                sounding_powerlaw_sens(hydro,wet)
            else
                sounding_powerlaw_sens_display
            end
        else
            % estimate the coefficients for individual interferograms
            sounding_powerlaw_sens(hydro,wet)
        end
 
        
    else
        fprintf('Based on specified sounding period \n\n')
        % When having sounding data estimate the powerlaw and b coefficient
        [alpha_log_all,alpha,h0,n_soundings] = sounding;


        % setting the estimated parameters in the parameter file
        fprintf('Updating the parm_aps list \n')
        setparm_aps('powerlaw_h0',h0)
        setparm_aps('powerlaw_alpha',alpha)
         
        % saving the data
        save([save_path filesep 'aps_p_step1.mat'],'alpha','h0','look_angle','lambda','sounding_h0','sounding_h_alpha_thres','sounding_start_date','sounding_end_date','sounding_time_stamp')

    end
end
if start_step==1 && strcmp(sounding_data,'n')
     fprintf('\n\n Step 1: powerlaw coefficents set by user \n\n')
     
     % as set by the user
     h0 = getparm_aps('powerlaw_h0');
     alpha = getparm_aps('powerlaw_alpha');
end

%% PART 2: Rotate the dataset to get mimumal white space around it
% Scaling and interpolation of the data to a regular grid
if start_step<=2 && end_step >=2
    fprintf('\n\nStep 2: Rotate the dataset \n')
    
    % getting the variables from the parm_aps file
    heading = getparm_aps('heading',1);
    crop_flag =  getparm_aps('crop_flag',1);

    % loading the lonlat information 
    ll_matfile = getparm_aps('ll_matfile',1);
    ll = load(ll_matfile);
    ll = ll.lonlat;
    
    % retrieving the coefficients
    h0 = getparm_aps('powerlaw_h0',1);
    alpha = getparm_aps('powerlaw_alpha',1);
    getparm_aps('powerlaw_xy_res',1);
    phuw_matfile = getparm_aps('phuw_matfile',1);
    hgt_matfile = getparm_aps('hgt_matfile',1);
    DEM_corr = getparm_aps('powerlaw_DEM_corr',1);
    bperp_matfile = getparm_aps('bperp_matfile',1);

    % rotating the data
    [xy_rotatelocal_full,xy_local_full,local_origin_full] = ll2rotatelocal(ll,heading);
    
    
    % cropping if needed
    crop = [];
    if strcmpi(crop_flag,'y')
        if exist('area_ex.mat','file')==2 
            crop = load('area_ex.mat'); 
        elseif exist('../area_ex.mat','file')==2
            crop = load('../area_ex.mat'); 
        else
            crop = [];
            fprintf(['There is no area_ex.mat file found, no data is cropped... \n'])    
        end
    end
    if ~isempty(crop)
        ix_remove =  inpolygon(ll(:,1),ll(:,2),crop.lonlat(:,1),crop.lonlat(:,2));
    else
        ix_remove = [];
    end
    ll(ix_remove,:) =[]; 
    
    % get the InSAR convexhull -- not needed at this stage!
    ix_convexhull = convhull(ll(:,1),ll(:,2));
    InSAR_convexhull = ll(ix_convexhull,:);
    if exist('tca_support.mat','file')==2
      save('tca_support.mat','-append','InSAR_convexhull')
    else
      save('tca_support.mat','InSAR_convexhull')        
    end
    clear InSAR_convexhull
    
    
    
    % rotating the data using all data including the region cropped. 
    % This will be used in the last stage of the power-law code.
    [xy_rotatelocal,xy_local,local_origin] = ll2rotatelocal(ll,heading);
    

    fprintf('        Scaling and interpolation to a regular grid \n\n')
    % loading the unwrapped interferograms
    ph = load(phuw_matfile);
    ph = ph.ph_uw; 
    
    % keep track of those pixels which are NaN
    ix_phase_nan = isnan(ph);
    
    % remove crop if needed
    ph(ix_remove,:)=[];
    
    
    % loading the topography information
    hgt = load(hgt_matfile);
    hgt = hgt.hgt; 
    % remove crop if needed
    hgt(ix_remove,:)=[];

    % geting the number of interferograms and points from the phase matrix
    n_interferograms= size(ph,2);
    n_points = size(ph,1);
    
    % Scaling of the topography according to the powerlaw
    % ifg based powerlaw correction?
    if length(alpha)==n_interferograms
        hgt_max = max(hgt);
        hgt_scaled=repmat(hgt,1,size(ph,2));
        ifg_based_correction='y';

        for k=1:length(alpha)
            if h0(k)*1000-hgt_max<=0
                error('myApp:argChk', ['Power law is not valid above h0 heights,... \nAbort,... \n'])       
            end
            hgt_scaled(:,k)=(h0(k)*1000-hgt_scaled(:,k)).^(alpha(k));
        end
    else
        % Scaling of the topography according to the powerlaw
        ifg_based_correction='n';
        ix = find(h0*1000-hgt < 0);
        if isempty(ix)~=1
            error('myApp:argChk', ['Power law is not valid above h0 heights,... \nAbort,... \n'])       
        end
        hgt_scaled = (h0*1000-hgt).^(alpha);
    end

    % define a new dataset were the height is actually the first column and
    % the interferograms are in the other columns
    % optional remove an estimate for the DEM errors 
    if strcmp(DEM_corr,'y') && n_interferograms>5;
       % these are erros that scale with perpendicular baseline
       % loading the perpendicualr baseline information
       bperp = load(bperp_matfile);
        if strcmp(stamps_processed,'y')
          bperp = bperp.bperp; 
       end
       % checking the size of bperp
       if size(bperp,2) > 1
          bperp = bperp';
          if size(bperp,2) >1
              error('myApp:argChk', ['bperp is not a vector,... \nAbort,... \n'])       
          end
       end
       % estimating the correlated errors
       DEM_corr_e = lscov(bperp,ph')';
       
       % removing DEM correlated errors 
       ph_temp = ph-repmat(bperp',n_points,1).*repmat(DEM_corr_e,1,n_interferograms);
       dataset = [hgt_scaled ph_temp];
       clear ph_temp A
    else
        if strcmp(DEM_corr,'y') && n_interferograms<=5
            fprintf('Not enough interferograms to make a reliable estimate for the DEM error \n')
            DEM_corr = 'n';
        end
        dataset = [hgt_scaled ph];
        DEM_corr_e = zeros([n_points 1]);
    end
    
    
    
    
    %% do an initial reference correction for each interferogram. 
    % This is to scope for regions that are small and which can have bandfiltering
    % arctifacts introduced. What happens is that each interferogram is
    % shifted such that the phase is zero at h0. Note that this is a long
    % wavelength shift that will be filtered out latter on. But in case of
    % irregularities it will reduce introduced arctifacts. To cope with
    % spatial varying signals this offset estimation is done over multiple
    % windows, but the constant offset is the same for all windows.     
    % get the position of some subwindows

    [iterate] = window_generation(xy_rotatelocal,1);
  
    counter=0;
    for k=1:length(iterate.window_ix)
        counter = counter + size(iterate.window_ix{k},1);
    end

    for k=1:n_interferograms
        % separate between ifg based correction or using a single alpha and h0
        if length(alpha)==n_interferograms
            height_data = dataset(:,k);
            scaling = 1./nanmean(abs(dataset(:,k)));
            phase_data = dataset(:,n_interferograms+k);
        else
            scaling = 1./nanmean(abs(dataset(:,1)));
            height_data = dataset(:,1);
            phase_data = dataset(:,1+k);
        end
        

        % estimating the reference, forced to be the same for all windows
        A_temp = zeros([counter length(iterate.window_ix)]);
        data_temp =  zeros([counter size(phase_data,2)]);
        counter =1;
        for kk=1:length(iterate.window_ix)
            A_temp(counter:counter+size(iterate.window_ix{kk},1)-1,kk) = [height_data(iterate.window_ix{kk},1)]*scaling;
            data_temp(counter:counter+size(iterate.window_ix{kk},1)-1,:) = [phase_data(iterate.window_ix{kk},1)];
            counter =  counter + size(iterate.window_ix{kk},1);

        end
        A_temp = [A_temp ones([size(data_temp,1) 1])];
        
        % allowing for NaN in interferogram {DB}
%        ix = find(isnan(data_temp(:,1))~=1);       
        ix = ~((isnan(data_temp(:,1)) + sum(isnan(A_temp),2))>=1);
        coeff = lscov(A_temp(ix,:),data_temp(ix,1));

        % plot the interferogram before and after the reference shift.
        if plot_flag ==1
           figure('name',['Correction of the reference for ifg ', num2str(k)])
           subplot(2,1,1)
           plot(height_data,phase_data,'k.')
           if length(alpha)==n_interferograms
                temp =getparm_aps('powerlaw_alpha');
                xlabel(['(h_0-h)^{' num2str(temp(k)) '}'])
           else
                xlabel(['(h_0-h)^{' num2str(getparm_aps('powerlaw_alpha')) '}'])              
           end
           ylabel('Phase')
        end

       if plot_flag ==1
           for kk=1:length(iterate.window_ix)
                y_temp = [[height_data(iterate.window_ix{kk},1)]*scaling ones(size([height_data(iterate.window_ix{kk},1)]*scaling))]*[coeff(kk) coeff(end)]';
                hold on
                plot([height_data(iterate.window_ix{kk},1)],y_temp,'r-')
           end       
           title('power-law plot with ifg reference to be the mean')
       end
       
       % correct the data for the reference
       data_temp = data_temp - coeff(end);
       % separate between ifg based correction or using a single alpha and h0
       if length(alpha)==n_interferograms
           dataset(:,n_interferograms+k) = dataset(:,n_interferograms+k) - coeff(end);
       else
           dataset(:,k+1) = dataset(:,k+1) - coeff(end);
       end
       
       % plot the corrected result
       if plot_flag ==1
           subplot(2,1,2)
           if length(alpha)==n_interferograms
               plot(dataset(:,k),dataset(:,n_interferograms+k),'k.')
               temp =getparm_aps('powerlaw_alpha');
               xlabel(['(h_0-h)^{' num2str(temp(k)) '}'])
           else
              plot(dataset(:,1),dataset(:,k+1),'k.')
              xlabel(['(h_0-h)^{' num2str(getparm_aps('powerlaw_alpha')) '}'])
           end
           ylabel('Phase shifted to h0 reference')
           
           % plot line on top
           coeff = lscov(A_temp,data_temp(:,1));
           for kk=1:length(iterate.window_ix)
                if length(alpha)==n_interferograms
                    y_temp = [[dataset(iterate.window_ix{kk},k)]*scaling ones(size([dataset(iterate.window_ix{kk},k)]*scaling))]*[coeff(kk) coeff(end)]';
                else
                    y_temp = [[dataset(iterate.window_ix{kk},1)]*scaling ones(size([dataset(iterate.window_ix{kk},1)]*scaling))]*[coeff(kk) coeff(end)]';
                end
                hold on
                plot([dataset(iterate.window_ix{kk},k)],y_temp,'r-')

           end
           title('power law plot with h0 as reference')
       end
    end

    
    %% Extrapolate such the convexhull contains the full grid
    xy_res = getparm_aps('powerlaw_xy_res');
    [xy_local_box,dataset_box] = extrapolate_local_new(xy_rotatelocal,dataset, xy_res(1), xy_res(2));

    % Interpolate to a regular grid
    [X_regular,Y_regular,Z_regular] = interpolate_regular(xy_local_box,dataset_box, xy_res(1), xy_res(2));
    
    % Saving the results of this step
    aps_save([save_path filesep 'aps_p_step2.mat'],xy_rotatelocal,xy_rotatelocal_full,ifg_based_correction,xy_local,xy_local_full,heading,ll,local_origin,local_origin_full,X_regular,Y_regular,Z_regular,xy_rotatelocal,xy_res,DEM_corr,DEM_corr_e,ix_remove,ix_phase_nan)
    
end


%% PART 3: Bandfiltering and interpolation back to a local grid
if start_step<=3 && end_step >=3
    fprintf('Step 3: Bandfiltering of the regular grid and converting back to a local grid \n\n')
    % Loading the data from the previous step
    if start_step==3
        load([save_path filesep 'aps_p_step2.mat'])
    end

    
     % loading the lonlat information 
    ll_matfile = getparm_aps('ll_matfile');
    temp = load(ll_matfile);
    if strcmp(stamps_processed,'y')
        temp = temp.lonlat;
    end
    
    n_points_original = size(temp,1);
    clear temp;
    
    % getting the extend of the spatial bandfilters
    spatial_bands = getparm_aps('powerlaw_spatial_bands');
    
    % Bandfiltering
    bandfiltering(Z_regular,xy_res(1),xy_res(2),spatial_bands,save_path,ifg_based_correction)
    
    % Interpolate back to local grid
    n_datasets = size(Z_regular,3);
    % calling the function for each bandfiltered dataset
    fprintf(['*Interpolate the band filtered data to the local grid: \n'])
    if strcmp(ifg_based_correction,'y')
        h_ifg_number = n_datasets/2;
    else
        h_ifg_number = 1;
    end
    for k=1:n_datasets
       
        % the file names, vary depending if the correction is varying for each interferogram       
        if k<=h_ifg_number && h_ifg_number~=1
            % thse are interferograms
            load_name = ['bandfilter_regular_hgt_ifg_' num2str(k) '.mat'];
            save_name = ['bandfilter_local_hgt_ifg' num2str(k) '.mat'];
        elseif k<=h_ifg_number && h_ifg_number==1
            % these are the heights
            load_name = 'bandfilter_regular_hgt.mat';
            save_name = 'bandfilter_local_hgt.mat';
        else
            % thse are interferograms
             load_name = ['bandfilter_regular_ifg_' num2str(k-h_ifg_number) '.mat'];
             save_name = ['bandfilter_local_ifg_' num2str(k-h_ifg_number) '.mat'];
        end

        dataset_band = load([save_path filesep load_name]);
        
        % Storing the information what type of filtering was done for which band
        dimension_filter = dataset_band.dimension_filter;      % 1 for 1D, 2 for 2D        
        % loading the data for each dataset
        dataset_band = dataset_band.data_band_out;
       
        % computing the bandfiltered data at the local grid
        % note that this is xy_rotatelocal and does not include the cropped
        % region!
        [z_local_out_temp] = interpolate2local(X_regular,Y_regular,dataset_band,xy_rotatelocal);
        % in case data has been cropped then path the cropped region with NaNs
        z_local_out = NaN([n_points_original size(z_local_out_temp,2)]);
        
        if isempty(ix_remove)
            z_local_out = z_local_out_temp;            
        else
            z_local_out(~ix_remove,:) = z_local_out_temp;
        end
        
        % Check in case some of the orginal phase values were nan.
        % if so set them to nan here too.
        if k<=h_ifg_number && h_ifg_number~=1
            % these are interferograms heights not the phase. No need to
            % modify this here.
        elseif k<=h_ifg_number && h_ifg_number==1
            % these are the heights, nothing needs to be done as the
            % heights remain the same for each interferogram.
        else
            % these are interferograms
            fprintf('Allowing for NaN''s in the interferograms \n')
            z_local_out(ix_phase_nan(:,k-h_ifg_number),:)=NaN;
        end
  
        
        % saving the data for each dataset
        aps_save([save_path filesep save_name],z_local_out,spatial_bands,dimension_filter);
        % outputting the progress to the user
        fprintf(['Progress: ' num2str(k) '/' num2str(n_datasets) ' done \n'])
    end
    
    % Saving the results of this step
    aps_save([save_path filesep 'aps_p_step3.mat'],ifg_based_correction,n_datasets,local_origin,local_origin_full,xy_rotatelocal,xy_rotatelocal_full,DEM_corr,DEM_corr_e,dimension_filter,ix_remove,ix_phase_nan)
end

%% PART 4: Estimating the scaling coefficient of the powerlaw for the local dataset 
% and computing the corresponding troposheric delay
if start_step<=4 && end_step >=4
    
    
    fprintf('\n\nStep 4: Estimating scaling coefficient of the powerlaw locally \n')
    fprintf('and estimating the tropospheric interferometric phase delay \n\n')
    
    % Loading the data from the previous step
    if start_step==4
        load([save_path filesep 'aps_p_step3.mat'])
    end

    % loading information from the parm file
    hgt_matfile = getparm_aps('hgt_matfile');
    powerlaw_all_bands = getparm_aps('powerlaw_all_bands');
    hgt = load(hgt_matfile);
    hgt = hgt.hgt; 
    
    % checking if its an interferogram based correction or not:
    if strcmp(ifg_based_correction,'y')
        n_interferograms = n_datasets/2;
    else
        n_interferograms = n_datasets-1;
    end
    
    
    % initialization
    n_points = size(hgt,1);
    ph_tropo_powerlaw = NaN([n_points n_interferograms]);
    K_tropo_powerlaw = NaN([n_points n_interferograms]);
    perc_cons_K_powerlaw= NaN([ n_interferograms 1]);
    if strcmp(powerlaw_all_bands,'y')
        ph_tropo_powerlaw_bands = NaN([n_points n_interferograms size(getparm_aps('powerlaw_spatial_bands'),1)]);        
        K_tropo_powerlaw_bands = NaN([n_points n_interferograms size(getparm_aps('powerlaw_spatial_bands'),1)]);
    end
    
    
    % Linear estimation for multiple windows based on bandfitlered data
    % final slope value is estimated for reliable bands only
    iterate = [];
    % getting the patch information from the parm_aps file
    
    
    % processing only those images that are not dropped in case of stamps
    ix_ifgs = 1:n_interferograms;
    if strcmp(stamps_processed,'y')
       ix_drop_ifgs = getparm('drop_ifg');
    else
        ix_drop_ifgs = [];
    end
    ix_ifgs(ix_drop_ifgs)=[];
    
    for k=1:length(ix_ifgs)
        
        ix_interferogram = ix_ifgs(k);

        % the file names, vary depending if the correction is varying for each interferogram       
        if strcmp(ifg_based_correction,'y') && n_interferograms~=1
            % the first dataset was set the be the topography
            topo_data = load([save_path filesep 'bandfilter_local_hgt_ifg' num2str(ix_interferogram) '.mat']);
            topo_data = topo_data.z_local_out;
            % the other dataset were the phase
            phase_data = load([save_path filesep 'bandfilter_local_ifg_' num2str(ix_interferogram) '.mat']);
            phase_data = phase_data.z_local_out;
            ifg_number=k;
        else
            % the first dataset was set the be the topography
            topo_data = load([save_path filesep 'bandfilter_local_hgt.mat']);
            topo_data = topo_data.z_local_out;
            % the other dataset were the phase
            phase_data = load([save_path filesep 'bandfilter_local_ifg_' num2str(ix_interferogram) '.mat']);
            phase_data = phase_data.z_local_out;
            ifg_number=1;
        end
                   
        % keep the figures of the outlier rejection as validation
        if k==0
            save_path_outlier_windows = [save_path filesep 'ifg_' num2str(ix_interferogram) '_outlier_window'];
        else
            save_path_outlier_windows=[];
        end
        
        if k==10000
            debug_fig_linear=1;
        else
            debug_fig_linear=0;
        end
        
        [ph_tropo_powerlaw_temp,ph_tropo_powerlaw_band_temp,iterate,K_temp,K_band,perc_cons_K] = aps_powerlaw_linear_local(xy_rotatelocal_full,local_origin_full,topo_data,phase_data,iterate,dimension_filter,save_path_outlier_windows,ifg_number,ix_phase_nan(:,ix_interferogram));
        
        
        % storing of the data
        % Computation of the tropospheric delay for all ifgs
        ph_tropo_powerlaw(:,ix_interferogram) = ph_tropo_powerlaw_temp;
        K_tropo_powerlaw(:,ix_interferogram) = K_temp;
        perc_cons_K_powerlaw(ix_interferogram,1) = perc_cons_K;
        if strcmp(powerlaw_all_bands,'y')
            % keeping all bands with in the columns the ifgs and the third
            % dimention the bands
            ph_tropo_powerlaw_bands(:,ix_interferogram,:) = ph_tropo_powerlaw_band_temp;
            K_tropo_powerlaw_bands(:,ix_interferogram,:)=K_band;
            
            % keeping the interferogram and the respective bands.
            eval(['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram) '= ph_tropo_powerlaw_band_temp;']);
            eval(['K_tropo_powerlaw_ifg_' num2str(ix_interferogram) '= K_band;']);
            
            clear ph_tropo_powerlaw_band_temp
            % saving the data from this step
            if strcmp(stamps_processed,'y')
                % This is StaMPS
                if strcmp(getparm('small_baseline_flag'),'y')
                    if exist(apsbandssbname,'file')==2
                        save(apsbandssbname,'-append',['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                    else
                        save(apsbandssbname,['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                    end
                else
                    if exist(apsbandsname,'file')==2
                        save(apsbandsname,'-append',['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                    else
                        save(apsbandsname,['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                    end
                end
            else
                % This is not StaMPS
                if exist(apsbandsname,'file')==2
                    save(apsbandsname,'-append',['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                else
                    save(apsbandsname,['ph_tropo_powerlaw_ifg_' num2str(ix_interferogram)],['K_tropo_powerlaw_ifg_' num2str(ix_interferogram)])
                end  
            end   
            eval(['clear ph_tropo_powerlaw_ifg_' num2str(ix_interferogram) ' K_tropo_powerlaw_ifg_' num2str(ix_interferogram)  ';']);
        end
        
        % give output to the screen
        fprintf(['Progress ifgs: ' num2str(k) '/' num2str(length(ix_ifgs)) ' done \n'])
    end
    % saving the data from this step
    aps_save([save_path filesep 'aps_p_step4.mat'],ph_tropo_powerlaw,DEM_corr,DEM_corr_e,iterate,K_tropo_powerlaw)
    
    % saving the data into the final variable togehter were other
    % correction values are saved.
    ph_tropo_powerlaw_or= ph_tropo_powerlaw;
    K_tropo_powerlaw_or=K_tropo_powerlaw;
    if strcmp(stamps_processed,'y')
        % This is StaMPS
        if strcmp(getparm('small_baseline_flag'),'y')
            if exist(apssbname,'file')==2
                save(apssbname,'-append','ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
            else
                save(apssbname,'ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
            end
            
            % saving the band information
            if strcmpi(powerlaw_all_bands,'y')
                if exist(apsbandssbname,'file')==2
                    save(apsbandssbname,'-append','ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
                else
                    save(apsbandssbname,'ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
                end
            end
        else
            % this is single master
            if exist(apsname,'file')==2
                save(apsname,'-append','ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
            else
                save(apsname,'ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
            end
            
            % saving the band information
            if strcmpi(powerlaw_all_bands,'y')
                if exist(apsbandsname,'file')==2
                    save(apsbandsname,'-append','ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
                else
                    save(apsbandsname,'ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
                end
            end
        end
    else
        % This is not StaMPS
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
        else
            save(apsname,'ph_tropo_powerlaw','K_tropo_powerlaw','perc_cons_K_powerlaw')
        end  
        
         % saving the band information
        if strcmpi(powerlaw_all_bands,'y')
            if exist(apsbandsname,'file')==2
                save(apsbandsname,'-append','ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
            else
                save(apsbandsname,'ph_tropo_powerlaw_or','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','K_tropo_powerlaw_bands','perc_cons_K_powerlaw')
            end
        end
    end
    setparm_aps('powerlaw_kept',[0]);

end

if start_step<=5 && end_step >=5 &&  exist('aps_RMSE_comparison')==2
    fprintf('\n\nStep 5: Powerlaw band comparison with reference\n')

    if strcmpi(getparm_aps('stamps_processed'),'y')
        % check which technique it needs to be compared with
        technique_number = -10;
        while isnumeric(technique_number)~=1 | (technique_number<1 | technique_number>7)
            technique_number = input(['*Select the technique number you want to compare with: \n---- MERIS (1), MODIS (2), MODIS Recal (3), ERA-I (4), WRF (5), u(sb) (6), u(sb)-d (7) : \n'],'s');
            technique_number = floor(str2num(technique_number));
        end
        % the technique component as a string
        if technique_number==1
            comp_str = 'a_mi';
        elseif technique_number==2
            comp_str = 'a_MI';
        elseif technique_number==3
            comp_str = 'a_RMI';
        elseif technique_number==4
            comp_str = 'a_e';
        elseif technique_number==5
            comp_str = 'a_w';
        elseif technique_number==6
            if strcmpi(getparm('small_baseline_flag'),'y')
                comp_str = 'usb';
            else
                comp_str = 'u';
            end
        elseif technique_number==7
            if strcmpi(getparm('small_baseline_flag'),'y')
                comp_str = 'usb-d';
            else
                comp_str = 'u-d';
            end
        end

        % check if we should add the hydrostatic delays too in case this is
        % possible
        if technique_number==1 | technique_number==2 | technique_number==3
            hydro_comp = inf;
            while isnumeric(hydro_comp)~=1 | (hydro_comp<0 | hydro_comp>2)
                hydro_comp = input(['*Do you want to add a hydrostatic component:\n---- None (0), ERA-I hydr (1), WRF hydro (2): \n'],'s');
                hydro_comp = floor(str2num(hydro_comp));
            end
        else
            hydro_comp = 0;
        end
        % the hydrostatic component as a string
        hydro_str = '';
        if hydro_comp==1
            hydro_str = '+a_eh';
        elseif hydro_comp==2
            hydro_str = '+a_wh';
        end

        % The reference technique is combination of both
        tech_ref = [comp_str hydro_str];
        clear hydro_str comp_str technique_number hydro_comp


        % running the RMSE script
        aps_RMSE_comparison(tech_ref,'a_pbands')
        str = '';
        while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
            str = input(['Do you want to update ''powerlaw_kept'' to a different band? [y/n] \n'],'s');
        end
        if strcmpi(str,'y')
            bandnumber = -10;
            while (bandnumber<0 | bandnumber>size(getparm_aps('powerlaw_spatial_bands'),1))
                bandnumber = input(['Which the band? (row number of the band to keep, or 0 for all bands combined) \n'],'s');
                bandnumber = str2num(bandnumber);
                if isempty(bandnumber)
                    bandnumber = -10;
                end
            end  
            setparm_aps('powerlaw_kept',[bandnumber]);
        end
    else
       fprintf('This option uses the StaMPS ps_plot function, and therefore only works if your data is stamps processed \n') 
    end

    
end
