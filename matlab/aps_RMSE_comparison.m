function [] = aps_RMSE_comparison(ref_technique,other_techniques,save_dir_path,save_dir_folder,units_flag,deramp_flag,fig_prop)
% Script that computes the RMSE between a reference technique and all the
% other technqiues. Optional is the profile flag, which when 1 will run the
% profiler code as well. Some hardoced values in the RMSE comparision code
% can still be changed in the code like, the percentage threshold of
% minimal common points, a deramping option, a reference option, and as
% last a option to crop out a specific region.
% Both the 'ref_technique' and 'other_techniques' variables are strings with 
% the abbreviated strings of the techniques. To specificy multiple
% techniques for the latter use tyhe & charakter, e.g. 'a_MI+a_eh&a_p'.
%
% The struct fig_prop can contain ylimits
%
% NOTE: the interferogram numbers displayed are the number of the
% interferograms kept. Not dropped. So in case of a dropped interferogram
% the number does not need to correspond. 
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
% By Bekaert David - University of Leeds.
% This approach is used in Bekaert et al., 2014 as comparison approach
% between techniques
%
%
% modifications
% 09/2014   DB:     Include the option tp call the profiling code directly
% 03/2015   DB:     change reference to black to difference more from the
%                   other methods. Change the plotting for powerlaw
%                   comparison of band.
%                   Make the interferogram plotting optional
% 03/2015   DB:     Include cm option as units. Plot the meris theoretical
%                   accuracy too when using MERIS as reference technique.
% 03/2015   DB:     Change the marker of the reference technqiue to a square
% 03/2015   DB:     Allow input of a matrix too
% 04/2015   DB:     Convert the theoretical accuracy of meris to LOS. 
%                   Difference between RMSE and RMS theoretical accuracy
% 04/2015   DB:     Add option to define the figure properties in advance.
% 02/2016   DB:     Include a fix when showing a single interferogram
%                   Fix for y-axis limits in case the correction is larger
%                   in RMSE. Fix in saving the figures.
% 08/2016   SSS:    compatibility with older matlab version for string splitting 


%% USER handles
% directory where the figures need to be saved together with the stats
save_dir = 'RMSE_comparison';
save_dir_bands = 'RMSE_bands';
save_dir_data = 'RMSE_own_data';

% radar wavelength in m converted to cm
lambda = getparm_aps('lambda')*100;     % cm units         
% retrieving the look angle
lookangle = getparm_aps('look_angle');
if ischar(lookangle) & exist(lookangle,'file')==2
    lookangle = load(lookangle);
    lookangle = lookangle.la;
end
lookangle = nanmean(lookangle);         % rad units


% meris uncertainty on PWV in cm units
PWV_meris = 1.1/10;                     % cm units
% compute the uncertancy of the MERIS delay for a difference between two
% epochs (reason for sqrt(2)). 
meris_uncertaincy_cm = 6.2*PWV_meris*sqrt(2)./cos(lookangle);
meris_uncertaincy_rad = meris_uncertaincy_cm*4*pi./lambda/cos(lookangle);


% threshold of the minimum percentage of points used to compute RMSE
perc_threshold = 50;


% when 1 the reference is set based on the mean of all available points in
% each ifgs.
reference_flag = 1;

% when 1 crop out a polygon. This needs to be defined as lonlat variable 
crop_file_path = '/nfs/a1/insar/mexico/envisat/track_255/SSE_model/programs/slip_model_code/model_cropout_MC.mat';
crop_out_flag = 0;



% the reference technique used to compare all the other data to
if nargin <1 | isempty(ref_technique)
    str = '';
    % compare with respect to the following techniques. 
    fprintf(['To indicate more use the & symbol.\n'])
    fprintf(['You can alsu use u option or a subtract flag of the u option.\n']);
    fprintf(['For the APS you need to indicate the type of correction using a_l, a_M+a_eh etc flag. \n']);
   
    while isempty(str) && isnumeric(str2num(str))==1
        str = input(['What will be the reference technique? (e.g. a_mi+a_eh) \n\n'],'s');
    end
    ref_technique = str;
    clear str
end
if nargin <2 | isempty(other_techniques)
    str = '';
    % compare with respect to the following techniques. 
    fprintf(['To indicate more use the & symbol.\n'])
    fprintf(['You can alsu use u option or a subtract flag of the u option.\n']);
    fprintf(['For the APS you need to indicate the type of correction using a_l, a_M+a_eh etc flag. \n']);
    while isempty(str) && isnumeric(str2num(str))==1
        str = input(['What will be the other technique(s)? (e.g. a_mi+a_eh) \n\n'],'s');
    end
    other_techniques = str;
    clear str
end


if nargin<3
    save_dir_path = [];
end
if nargin<4
    save_dir_folder = [];
end
if nargin<5
    units_flag='rad';
end

if ~(strcmpi(units_flag,'rad') || strcmpi(units_flag,'cm'))
    error('Give a valid unit flag, either a string ''rad'' or ''cm''')
end

if nargin<6
    % when 1 deramp all images else do not
    deramp_flag = 0;    
end

% related to figure properties
if nargin <7
    fig_prop=[];
end
if isfield(fig_prop,'ylimits')
    ylimits = fig_prop.ylimits;
else
    ylimits=[];
end
if isfield(fig_prop,'current_colors')
    current_colors = fig_prop.current_colors;
else
    current_colors=[];
end

    

%% actual code

% paths
InSAR_datapath = pwd;
curdir = pwd;


% getting the original reference area 
ref_radius_original = getparm('ref_radius');
setparm('ref_radius',nan)


% Checking if this is to plot the powerlaw bands or actual different
% correction methods
compare_bands=0;
if ischar(other_techniques)
    if strcmpi(other_techniques,'a_pbands')
        compare_bands = 1;
        save_dir= save_dir_bands;
    else
        compare_bands = 0;
    end
else
    % this is based on a matrix specified
    other_techniques_data = other_techniques;
    clear other_techniques
    other_techniques = repmat(' &',1,size(other_techniques_data,3));
    other_techniques=other_techniques(1:end-1);
    save_dir = save_dir_data;
    
    % set reference area of the data
    ix_ref = ps_setref;
    for k_temp = 1:size(other_techniques_data,3)
        other_techniques_data(:,:,k_temp) = other_techniques_data(:,:,k_temp)-repmat(nanmean(other_techniques_data(ix_ref,:,k_temp),1),size(other_techniques_data,1),1);
    end
end
% check if the reference technique is given as a string
if ~ischar(ref_technique)
    save_dir = save_dir_data;
    ref_technique_data = ref_technique;
    clear ref_technique;
    ref_technique = ' ';        % this means user has inputed his own dataset
    
    % set reference area of the data
    ix_ref = ps_setref;
    ref_technique_data = ref_technique_data-repmat(nanmean(ref_technique_data(ix_ref,:),1),size(ref_technique_data,1),1);
end    


% figure options
fontsize = 15;
if ~isempty(save_dir_folder)
   save_dir = save_dir_folder; 
end
if isempty(save_dir_path)
    save_path = [InSAR_datapath filesep save_dir];
else
    save_path = [save_dir_path filesep save_dir];
end
if exist(save_path,'dir')~=7
    mkdir(save_path)
else
    str = '';
    while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
        str = input(['This folder already exist, want to clear it? [y/n] \n'],'s');
    end
    if strcmp(str,'y')
        delete([save_path filesep '*'])
    end
end


% get if its small baseline or not 
small_baseline_flag = getparm('small_baseline_flag');
load('psver');
ps = load(['ps' num2str(psver) '.mat']);
if strcmp(small_baseline_flag,'y')
    aps_str = 'asb';
    tca_name = 'tca_sb2.mat';
    tca_bands_name = 'tca_bands_sb2.mat';
    
    % % for baselines you need ps of PS directory 
    if exist(['..' filesep 'ps1.mat'],'file')~=2
       bperp = NaN([length(ps.ifgday_ix) 2]) ; 
    else
        bperp = load(['..' filesep 'ps1.mat'],'bperp');
        bperp = bperp.bperp;
        bperp = bperp(ps.ifgday_ix);
    end
    dates=ps.ifgday;   

    if ~isempty(strfind(ref_technique,'u')) & isempty(strfind(ref_technique,'usb'))
       fprintf('You are displaying a phase type for PS not SB!\nThis has been changed to SB instead!\n') 
       ref_technique = ['usb' ref_technique(2:end)];
       fprintf([ref_technique '\n\n'])
    end
    if ~isempty(strfind(other_techniques,'u')) & isempty(strfind(other_techniques,'usb'))
       fprintf('You are displaying a phase type for PS not SB!\nThis has been changed to SB instead!\n') 
       ix = strfind(other_techniques,'u');
       other_techniques = [other_techniques(1:ix) 'sb' other_techniques(ix+1:end)];
       fprintf([other_techniques '\n\n'])
    end    
    ix_ifgs_keep = 1:ps.n_ifg;
else
    aps_str = 'a';   
    tca_name = 'tca2.mat';
    tca_bands_name = 'tca_bands2.mat';

    dates = [ps.day   repmat(ps.master_day,length(ps.day),1)]; 
    bperp = [ps.bperp zeros([length(ps.day) 1])];
    ix_ifgs_keep = 1:ps.n_ifg;
end

% remove that that have not been unwrapped before
ix_dropped = getparm('drop_ifg');
if strcmp(small_baseline_flag,'n')
    ix_dropped = unique([ix_dropped ps.master_ix]);
end
ix_ifgs_keep(ix_dropped)=[];
bperp(ix_dropped,:)=[];
dates(ix_dropped,:)=[];
    
% keep the original network for reference
dates_all = dates;
bperp_all = bperp;
    

% technique to be compared to
try 
    other_techniques = strsplit(other_techniques,'&');  %SSS 2/8/16
catch
    other_techniques=strread(other_techniques,'%s','delimiter','&')  %alternative to "strsplit" which is not availabe on matlab versions <2013a. SSS 2/8/16
end

% other techniques. This could also be bandfiltered powerlaw data
if compare_bands==1
    bandfilter_dataset = getparm_aps('powerlaw_spatial_bands');
    n_techniques = size(bandfilter_dataset,1);
else
    n_techniques = length(other_techniques);
end
legend_str = [];



if compare_bands==1
    n_ifgs = length(ix_ifgs_keep);
   

    % the legend strings are the actual band filters
    lenged_str = cell(n_techniques+1,1);
    for technique_counter=1:n_techniques
        if technique_counter==1
            legend_str{technique_counter+1}= (['Power-law (w+h) band: ' num2str([bandfilter_dataset(technique_counter,1)./1000]) ' - ' num2str([bandfilter_dataset(technique_counter,2)./1000]) ' km']);            
        else
            legend_str{technique_counter+1}= ([num2str([bandfilter_dataset(technique_counter,1)./1000]) ' - ' num2str([bandfilter_dataset(technique_counter,2)./1000]) ' km']);
        end
        save_str{technique_counter+1} = (['Powerlaw_WH_band_' num2str([bandfilter_dataset(technique_counter,1)./1000]) '_' num2str([bandfilter_dataset(technique_counter,2)./1000]) 'km']);

    end
    
    
    % initializing the data
    data = NaN([ps.n_ps n_ifgs n_techniques+1]);
    data_temp = load(tca_bands_name);
    data(:,:,2:end) = data_temp.ph_tropo_powerlaw_bands(:,ix_ifgs_keep,:);
    clear data_temp temp_keep

    % check for NaNs, normally all bands should be fine
    for technique_counter=1:n_techniques
        ix_no_correction_temp = find((sum(isnan(data(:,:,technique_counter+1))==0,1)==0)==1);
        data(:,ix_no_correction_temp,technique_counter+1) = NaN;
        ix_no_correction{technique_counter+1} = ix_no_correction_temp;
        clear ix_no_correction_temp        
    end
    
    % loading the data of the reference
    if ref_technique(1)=='a'   
        ps_plot([aps_str],ref_technique,-1,0,0,ix_ifgs_keep);

        % loading the data
        data_temp = load(['ps_plot_' aps_str  '.mat']);
        delete(['ps_plot_' aps_str '.mat']);

        % check if there is actually a delay otherwize make it all NaN
        ix_no_correction_temp = find((sum(data_temp.ph_disp~=0,1)==0)==1);
        data_temp.ph_disp(:,ix_no_correction_temp) = NaN;
        ix_no_correction_temp = find((sum(isnan(data_temp.ph_disp)==0,1)==0)==1);
        ix_no_correction{1} = ix_no_correction_temp;
        clear ix_no_correction_temp
        data(:,:,1) = data_temp.ph_disp;

        % getting the legend string
        [legend_str{1} save_str{1}]= aps_name(ref_technique);

        
                
    elseif ref_technique(1)=='u'
        ps_plot([ref_technique],-1,0,0,ix_ifgs_keep);
        % loading the data
        data_temp = load(['ps_plot_' ref_technique '.mat']);
        delete(['ps_plot_' ref_technique '.mat']);

        % check if there is actually a delay otherwize make it all NaN
        ix_no_correction_temp = find((sum(data_temp.ph_disp~=0,1)==0)==1);
        data_temp.ph_disp(:,ix_no_correction_temp) = NaN;
        ix_no_correction_temp = find((sum(isnan(data_temp.ph_disp)==0,1)==0)==1);
        ix_no_correction{1} = ix_no_correction_temp;
        clear ix_no_correction_temp
        data(:,:,1) = data_temp.ph_disp;

        
        % legend string
        legend_str{1} = ref_technique;
        save_str_temp = ref_technique;
        save_str_temp(save_str_temp=='-')='_';
        save_str{1}=save_str_temp;
        clear save_str_temp
        
    elseif ref_technique(1)==' '
        % check if there is actually a delay otherwize make it all NaN
        ix_no_correction_temp = find((sum(ref_technique_data~=0,1)==0)==1);
        ref_technique_data(:,ix_no_correction_temp) = NaN;
        ix_no_correction_temp = find((sum(isnan(ref_technique_data)==0,1)==0)==1);
        ix_no_correction{1} = ix_no_correction_temp;
        clear ix_no_correction_temp
        data(:,:,1) = ref_technique_data;
        
        % legend string
        legend_str{1} = 'own-data';
        save_str_temp = 'own-data';
        save_str_temp(save_str_temp=='-')='_';
        save_str{1}=save_str_temp;
        clear save_str_temp
        
    end  
    n_datasets = n_techniques+1;

    
    
else
    for technique_counter =1:n_techniques+1
        if technique_counter==n_techniques+1
            technique_str = ref_technique;
        else
           technique_str = other_techniques{technique_counter};
        end

        % Telling user which technique is being done
        fprintf(['Processing: ' technique_str '\n' ])

        % plotting the data
        % this is APS
        if technique_str(1)=='a'   
            ps_plot([aps_str],technique_str,-1,0,0,ix_ifgs_keep);

            % loading the data
            data_temp = load(['ps_plot_' aps_str '.mat']);
            delete(['ps_plot_' aps_str '.mat']);

            % check if there is actually a delay otherwize make it all NaN
            ix_no_correction_temp = find((sum(data_temp.ph_disp~=0,1)==0)==1);
            data_temp.ph_disp(:,ix_no_correction_temp) = NaN;

            ix_no_correction_temp = find((sum(isnan(data_temp.ph_disp)==0,1)==0)==1);


            % legend string
            [legend_str_temp save_str_temp]= aps_name(technique_str);

        % this is a phase type
        elseif technique_str(1)=='u'
            ps_plot([technique_str],-1,0,0,ix_ifgs_keep);
            % loading the data
            data_temp = load(['ps_plot_' technique_str '.mat']);
            delete(['ps_plot_' technique_str '.mat']);

            % check if there is actually a delay otherwize make it all NaN
            ix_no_correction_temp = find((sum(data_temp.ph_disp~=0,1)==0)==1);
            data_temp.ph_disp(:,ix_no_correction_temp) = NaN;
            ix_no_correction_temp = find((sum(isnan(data_temp.ph_disp)==0,1)==0)==1);


            % legend string
            legend_str_temp = technique_str;
            save_str_temp = technique_str;
            save_str_temp(save_str_temp=='-')='_';
            
        % this is when the user gave its own data matrix
        elseif  technique_str(1)==' '
            if technique_counter==n_techniques+1
                data_temp.ph_disp = ref_technique_data;
            else
                data_temp.ph_disp = other_techniques_data(:,:,technique_counter);
            end
            
            % check if there is actually a delay otherwize make it all NaN
            ix_no_correction_temp = find((sum(data_temp.ph_disp~=0,1)==0)==1);
            data_temp.ph_disp(:,ix_no_correction_temp) = NaN;
            ix_no_correction_temp = find((sum(isnan(data_temp.ph_disp)==0,1)==0)==1);
            
            % legend string
            legend_str_temp = 'own-data';
            save_str_temp = 'own-data';
            save_str_temp(save_str_temp=='-')='_';
            
        end

        % storing the data   
        if technique_counter==1
            data = NaN([size(data_temp.ph_disp)  n_techniques+1]);
        end
        % putting the reference data in the first layer
        if technique_counter==n_techniques+1
            data(:,:,1) = data_temp.ph_disp;
            legend_str{1} = legend_str_temp;
            save_str{1} = save_str_temp;
            ix_no_correction{1} = ix_no_correction_temp;
        % putting the other data after the reference, keep the same order
        % as given by the user
        else
            data(:,:,technique_counter+1) = data_temp.ph_disp;
            legend_str{technique_counter+1} = legend_str_temp;
            save_str{technique_counter+1} = save_str_temp;
            ix_no_correction{technique_counter+1} = ix_no_correction_temp;

        end
        clear data_temp ix_no_correction_temp
    end
    n_datasets = n_techniques+1;
end
n_ifgs = size(data,2);
setparm('ref_radius',ref_radius_original);

% crop out a region, by filling them with NaNs
if crop_out_flag==1
    poly = load(crop_file_path);
    ix_crop_out = inpolygon(ps.lonlat(:,1),ps.lonlat(:,2),poly.lonlat(:,1),poly.lonlat(:,2));
    for dataset_counter=1:n_datasets;
        data(ix_crop_out,:,dataset_counter) = NaN;
    end
    
    clear ix_crop_out
    
else
    crop_file_path = [];
end


% de-ramp the data when required
if deramp_flag==1
    for dataset_counter=1:n_datasets;
        data(:,:,dataset_counter) = ps_deramp(ps,data(:,:,dataset_counter));
    end
end


% convert the units when requested:
if strcmpi(units_flag,'cm')
    data = data.*lambda./4./pi;
end

% comparing datasets - two options implemented
% option 1: all datasets have the same points for a specific interferogram
% option 2: each dataset only compares points that it has mutual with the reference

% implementation of option 1 
data_option1 = data;
for dataset_counter=1:n_datasets;
    for ifg_counter=1:n_ifgs
        % take all the datasets for a specific interferogram
        temp = squeeze(data_option1(:,ifg_counter,:));
        % find if there are pixels not present in some datasets. If so put
        % all of them to a nan
        ix_nan = sum(isnan(temp),2)>=1;
        
        % putting NaN's everywhere for the interferogram for all datasets
        data_option1(ix_nan,ifg_counter,:)=NaN;
         
        % check if the reference needs to be computed based on all remaining points
        % in the interferogram
        if reference_flag==1
            % removing the mean for each interfergram, but all datasets same time
            temp =  squeeze(data_option1(:,ifg_counter,:));
            temp = temp- repmat(nanmean(temp,1),size(temp,1),1);
            % storing it back in the orginal data matrix
            data_option1(:,ifg_counter,:)=temp;
            clear temp
        end
    end      
end

% Computing the RMSE
residual_sq_option1 = (data_option1 - repmat(data_option1(:,:,1),[1 1 n_datasets])).^2;
RMSE_option1 =  sqrt(nanmean((residual_sq_option1),1));
% replace the first RMSE with the RMSS (Root Mean Squared Signal) 
RMSE_option1(1,:,1) = sqrt(nanmean(data_option1(:,:,1).^2,1));
% RMSE with first column the reference
RMSE_option1 = squeeze(RMSE_option1);           % (colum vector)
if n_ifgs==1 & size(RMSE_option1,1)>1
    RMSE_option1 = RMSE_option1';
end
% fraction of data used to do the RMSE computation
RMSE_option1_perc=round((sum(isnan(residual_sq_option1(:,:,1))~=1,1)./size(residual_sq_option1(:,:,1),1)).*100);

% checking if the RMSE is forfulling the percentage threshold of points
% used in the calculation
ix_failed =  RMSE_option1_perc<=perc_threshold;
ix_no_correction_temp = [1:n_ifgs];
ix_no_correction_temp=ix_no_correction_temp(ix_failed);
RMSE_option1(ix_failed,:)=NaN; 
temp_failed = setxor(ix_no_correction_temp,unique(cell2mat(ix_no_correction)));
if length(temp_failed)>0
    fprintf([num2str(length(temp_failed)) ' interferograms failed the percentage threshold \n'])
end
ix_no_correction{n_techniques+2}=ix_no_correction_temp;
clear ix_no_correction_temp



% getting the interferograms where RMSE where computed
ifgs_good = isnan(sum(RMSE_option1,2))~=1;
ix_ifgs = 1:n_ifgs;
ix_ifgs = ix_ifgs(ifgs_good);


% original ifgs list
temp = 1:ps.n_ifg;
temp = temp(ix_ifgs_keep);
fprintf(['original interferogram list: ' num2str(temp(ix_ifgs))  '\n'])


% setting the colorbar of the plots
if isempty(current_colors)
    h_test = figure;
    colormap default
    % getting the current colors and keep that for the next plot
    current_colors = colormap;
    ix_current_colors = round([1:n_techniques+1]'.*size(current_colors,1)./(n_techniques+2));
    current_colors = current_colors(ix_current_colors,:);
    current_colors(1,:)=[0 0 0];    % put the reference technique to black
    close(h_test)
end

% plotting the bar chart for each technique and interferogram
h_bar = figure('name','RMSE bar chart','position',[-6         440        1171         367]);
colormap(current_colors);
if n_ifgs==1
    dummy = [RMSE_option1(ix_ifgs,:); NaN.*ones(size(RMSE_option1(ix_ifgs,:)))];
    bar(dummy);
else
    bar(RMSE_option1(ix_ifgs,:));
end
legend(legend_str,'location','northoutside')
set(gca,'xtick',[1:sum(ifgs_good)],'xticklabel',num2str(ix_ifgs'))
xlabel('Interferogram number','fontsize',fontsize)
% check the legend based on the units selected
if strcmpi(units_flag,'cm')
    ylabel('RMSE [cm]','fontsize',fontsize)
else
    ylabel('RMSE [rad]','fontsize',fontsize)
end
set(gca,'fontsize',fontsize)
% take the maximum of the 5 techniques for each interferogram and add a
% small shift up for the text, but keep the shift constant.
if isempty(ix_ifgs) || isnan(nanmean(nanmean(RMSE_option1)))
    y_text = 1;
else
    y_text = max(RMSE_option1(ix_ifgs,:),[],2) + mean(max(RMSE_option1(ix_ifgs,:),[],2))*0.25;
end
ylim([0 max(y_text)+0.25*max(y_text)])
for ifgs_counter = 1:sum(ifgs_good)
    text(ifgs_counter-0.25, y_text(ifgs_counter), [ num2str(RMSE_option1_perc(ix_ifgs(ifgs_counter))) '%' ], 'VerticalAlignment', 'top', 'FontSize', fontsize)
end
RMSE_perc = RMSE_option1_perc(ix_ifgs);
set(h_bar,'PaperPositionMode','auto')
print(h_bar,'-depsc',[save_path filesep 'RMSE_bar_chart.eps'])
print(h_bar,'-dpng',[save_path filesep 'RMSE_bar_chart.png'])


% plotting for each technique all the interferogram on a line
scatterplot = [];           % saving the data for replotting
scatterplot_mean = [];  
% for the band comparison also include band information on axis
if compare_bands==1
    h_scatter = figure('name','RMSE scatter plot of all ifgs','position',  [ 127   257   797   576]);
    % the bands
    h1 = axes('position',[0.22 0.4 0.65 0.45]);    
    for technique_counter=1:n_techniques
        for ifg_counter=1:length(ix_ifgs)
            plot(mean(bandfilter_dataset(technique_counter,:)./1000),[RMSE_option1(ix_ifgs(ifg_counter),technique_counter+1)],'k^','markerfacecolor',current_colors(technique_counter+1,:))
            hold on
        end
        hold on
        plot(mean(bandfilter_dataset(technique_counter,:)./1000),[nanmedian(RMSE_option1(ix_ifgs,technique_counter+1))],'k^','markerfacecolor',current_colors(technique_counter+1,:),'markersize',15,'linewidth',2)

        % store data for saving
        scatterplot = [scatterplot ; [technique_counter.*ones([sum(ifgs_good) 1]) RMSE_option1(ix_ifgs,technique_counter+1)]];
        scatterplot_mean = [scatterplot_mean; [technique_counter,nanmedian(RMSE_option1(ix_ifgs,technique_counter+1))]];
    end
    % plotting the separation between RMS signal and RMS error
    hold on 
    set(gca,'YAxisLocation','right')
    xlim([0 ceil(max(max(bandfilter_dataset))./1000)])
    % in case MERIS is the reference plot the theoretical accuracy as well.
    if length(ref_technique)>2 && strcmp(ref_technique(1:3),'a_m')
        xlimits_extend = get(gca,'xlim');
        if strcmpi(units_flag,'cm')
            % technique error is about zero
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend)),'k-')
            % technique error is same as meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend))*sqrt(2),'k--')
            % technique error is twice meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend))*sqrt(5),'k:')
        else
            % technique error is about zero
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend)),'k-')
            % technique error is same as meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend))*sqrt(2),'k--')
            % technique error is twice meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend))*sqrt(5),'k:')
        end
    end
    % check the legend based on the units selected
    if strcmpi(units_flag,'cm')
        ylabel('RMS [cm]','fontsize',fontsize)
    else
        ylabel('RMS [rad]','fontsize',fontsize)
    end
    % setting the other part of the axis.
    if isempty(ylimits)
        ylimits = [0 ceil(nanmax(nanmax(RMSE_option1)))];
    end
    ylim(ylimits)
    set(gca,'fontsize',fontsize)
    set(gca,'XAxisLocation','top')     
    xlabel('Spatial wavelength [km]','fontsize',fontsize)
    ylabel('RMSE [rad]','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    yticks = get(gca,'ytick');
    ylimits = get(gca,'ylim');
   
    
    % plotting the reference
    h1 = axes('position',[0.1 0.4 0.1 0.45]);
    plot(repmat(-10,length(ix_ifgs),1),[(RMSE_option1(ix_ifgs,1))],'ksq','markerfacecolor',current_colors(1,:))
    hold on
    plot(-10,[nanmedian(RMSE_option1(ix_ifgs,1))],'ksq','markerfacecolor',current_colors(1,:),'markersize',15,'linewidth',2)
    scatterplot_ref = [repmat(-10,length(ix_ifgs),1),[(RMSE_option1(ix_ifgs,1))]];
    scatterplot_ref_mean = [-10 ; nanmedian(RMSE_option1(ix_ifgs,1))];
    % in case MERIS is the reference plot the theoretical accuracy as well.
    if length(ref_technique)>2 && strcmp(ref_technique(1:3),'a_m')
        xlimits_extend = get(gca,'xlim');
        if strcmpi(units_flag,'cm')
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend)),'k--','linewidth',2)
        else
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend)),'k--','linewidth',2)
        end
    end
    % check the legend based on the units selected
    if strcmpi(units_flag,'cm')
        ylabel('RMS [cm]','fontsize',fontsize)
    else
        ylabel('RMS [rad]','fontsize',fontsize)
    end
    set(h1,'ytick',yticks)
    ylim(ylimits)
    xticks = [-10];
    xticks_label = 'REF';
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',xticks_label);
    set(gca,'fontsize',fontsize)

    % plotting the bandwidths
    h2 = axes('position',[0.22 0.175 0.65 0.2]); 
    [sorted_mean, ix] = sort(nanmedian(RMSE_option1(:,2:end),1));
    temp = reshape(RMSE_option1,[],1);
    temp(isnan(temp))=[];
    temp = sort(temp,'descend');
    for counter=1:n_techniques
            technique_counter = ix(counter);
            plot(bandfilter_dataset(technique_counter,:)./1000,    counter.*[1 1] ,'k-','linewidth',2,'color',current_colors(technique_counter+1,:))
            hold on
    end
    xlim([0 ceil(max(max(bandfilter_dataset))./1000)])
    ylim([0 counter+1])
    set(gca,'YAxisLocation','right')
    xlabel('Spatial wavelength [km]','fontsize',fontsize)
    set(gca,'ytick',[]);
    set(gca,'yticklabel','');
    ylabel({'<- smaller','RMSE'},'fontsize',fontsize)
    set(gca,'fontsize',fontsize);
    hold on
    box on
    
    % saving the figure
    set(h_scatter,'PaperPositionMode','auto')
    print(h_scatter,'-depsc',[save_path filesep 'RMSE_scatter_plot.eps'])
    print(h_scatter,'-dpng',[save_path filesep 'RMSE_scatter_plot.png'])
    
else
    % technique comparison
    h_scatter = figure('name','RMSE scatter plot of all ifgs','position',[ 127   257   797   576]);
    h1 = axes('position',[0.22 0.1 0.65 0.45]);    
    for technique_counter=1:n_techniques
        plot(technique_counter,nanmedian(RMSE_option1(ix_ifgs,technique_counter+1)),'k^','markerfacecolor',current_colors(technique_counter+1,:),'markersize',15,'linewidth',2)
        hold on
    end
%     legend(legend_str{2:end},'location','northoutside')
    for technique_counter=1:n_techniques
        hold on
        plot(technique_counter.*ones([sum(ifgs_good) 1]) ,RMSE_option1(ix_ifgs,technique_counter+1)','k^','markerfacecolor',current_colors(technique_counter+1,:))
        hold on
        plot(technique_counter,nanmedian(RMSE_option1(ix_ifgs,technique_counter+1)),'k^','markersize',15,'linewidth',2)

        % store data for saving
        scatterplot = [scatterplot ; [technique_counter.*ones([sum(ifgs_good) 1]) RMSE_option1(ix_ifgs,technique_counter+1)]];
        scatterplot_mean = [scatterplot_mean; [technique_counter,nanmedian(RMSE_option1(ix_ifgs,technique_counter+1))]];
    end
    xlim([0 n_techniques+1])
    % in case MERIS is the reference plot the theoretical accuracy as well.
    if length(ref_technique)>2 && strcmp(ref_technique(1:3),'a_m')
        xlimits_extend = get(gca,'xlim');
        if strcmpi(units_flag,'cm')
            % technique error is about zero
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend)),'k-')
            % technique error is same as meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend))*sqrt(2),'k--')
            % technique error is twice meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend))*sqrt(5),'k:')
        else
            % technique error is about zero
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend)),'k-')
            % technique error is same as meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend))*sqrt(2),'k--')
            % technique error is twice meris error
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend))*sqrt(5),'k:')
        end
    end
    set(gca,'xtick',[],'xticklabel','')
    xlabel('Techniques','fontsize',fontsize)
    set(gca,'YAxisLocation','right')
    if isempty(ylimits)
        ylimits = [0 ceil(max(get(gca,'ylim')))];
    end
    ylim(ylimits)
    set(gca,'fontsize',fontsize)
    % check the legend based on the units selected
    if strcmpi(units_flag,'cm')
        ylabel('RMSE [cm]','fontsize',fontsize)
    else
        ylabel('RMSE [rad]','fontsize',fontsize)
    end
    yticks = get(gca,'ytick');
    
    
    % plotting the reference
    h1 = axes('position',[0.1 0.1 0.1 0.45]);
    plot(repmat(-10,length(ix_ifgs),1),[(RMSE_option1(ix_ifgs,1))],'ksq','markerfacecolor',current_colors(1,:))
    hold on
    plot(-10,[nanmedian(RMSE_option1(ix_ifgs,1))],'ksq','markerfacecolor',current_colors(1,:),'markersize',15,'linewidth',2)
    scatterplot_ref = [repmat(-10,length(ix_ifgs),1),[(RMSE_option1(ix_ifgs,1))]];
    scatterplot_ref_mean = [-10 ; nanmedian(RMSE_option1(ix_ifgs,1))];

    
    % in case MERIS is the reference plot the theoretical accuracy as well.
    if length(ref_technique)>2 && strcmp(ref_technique(1:3),'a_m')
        xlimits_extend = get(gca,'xlim');
        if strcmpi(units_flag,'cm')
            hold on
            plot(xlimits_extend,meris_uncertaincy_cm*ones(size(xlimits_extend)),'k--','linewidth',2)
        else
            hold on
            plot(xlimits_extend,meris_uncertaincy_rad*ones(size(xlimits_extend)),'k--','linewidth',2)
        end
    end
    % check the legend based on the units selected
    if strcmpi(units_flag,'cm')
        ylabel('RMS [cm]','fontsize',fontsize)
    else
        ylabel('RMS [rad]','fontsize',fontsize)
    end
    set(h1,'ytick',yticks)
    ylim(ylimits)
    xticks = [-10];
    xticks_label = 'REF';
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',xticks_label);
    set(gca,'fontsize',fontsize)

    
    % plotting the legend as two parts:
    n_legend = ceil((n_techniques+1)/2);
    % part 1
    h3 = axes('position',[0.25 0.55 0.3 0.3]);    
    for technique_counter=1:n_legend 
        % reference has a different marker
        if technique_counter==1
            plot(technique_counter,-100,'ksq','markerfacecolor',current_colors(technique_counter,:),'markersize',15,'linewidth',2)
        else
            plot(technique_counter,-100,'k^','markerfacecolor',current_colors(technique_counter,:),'markersize',15,'linewidth',2)
        end
        hold on
    end
    xlim([-1 1])
    ylim(ylimits)
    h2_legend = legend(legend_str{1:n_legend});
    set(h2_legend,'location','east')
    set(gca,'visible','off')
    set(gca,'fontsize',fontsize);
    
    % part 2
    h2 = axes('position',[0.55 0.55 0.3 0.3]);    
    for technique_counter=n_legend+1:n_techniques+1 
        plot(technique_counter,-100,'k^','markerfacecolor',current_colors(technique_counter,:),'markersize',15,'linewidth',2)
        hold on
    end
    xlim([-1 1])
    ylim(ylimits)
    h3_legend= legend(legend_str{n_legend+1:n_techniques+1});
    set(h3_legend,'location','west')
    set(gca,'visible','off')
    set(gca,'fontsize',fontsize);

    % saving the figure
    set(h_scatter,'PaperPositionMode','auto')
    print(h_scatter,'-depsc',[save_path filesep 'RMSE_scatter_plot.eps'])
    print(h_scatter,'-dpng',[save_path filesep 'RMSE_scatter_plot.png'])
end

if sum(sum(isnan(bperp_all)))==0
    % plotting the network of the interferograms used in the RMSE computation
    h_baselineplot = figure('name','Processed network');
    % plotting all ifgs that were considered
    for ifgs_counter=1:size(bperp_all,1)
        hold on
        plot([dates_all(ifgs_counter,1)  dates_all(ifgs_counter,2)],  [bperp_all(ifgs_counter,1) bperp_all(ifgs_counter,2)] ,'k-','linewidth',1) 
    end
    % plotting the network for which we have an APS correction 
    for ifgs_counter=1:sum(ifgs_good)
        hold on
        plot([dates(ix_ifgs(ifgs_counter),1)  dates(ix_ifgs(ifgs_counter),2)],  [bperp(ix_ifgs(ifgs_counter),1) bperp(ix_ifgs(ifgs_counter),2)] ,'k-','linewidth',2) 


        text(mean(dates(ix_ifgs(ifgs_counter),:)),mean(bperp(ix_ifgs(ifgs_counter),:))+20,num2str(ix_ifgs(ifgs_counter)),'fontsize',fontsize-2,'backgroundColor',[1 1 1],'Margin',0.01)
    end
    hold on
    % [dates_unique,ix_unique] = unique(dates(ix_ifgs,:));
    % bperp_unique = bperp(ix_ifgs,:);
    % bperp_unique = bperp_unique(ix_unique);
    [dates_unique,ix_unique] = unique(dates_all);
    bperp_unique = bperp_all;
    bperp_unique = bperp_unique(ix_unique);
    plot(dates_unique,bperp_unique,'ko','markerfacecolor','r','markersize',7)
    clear dates_unique ix_unique bperp_unique
    hold on
    % set(gca,'XTick',dates_num)
    datetick('x','mmm yy')
    set(gca,'fontsize',fontsize)
    box on
    ylabel('Bperp [m]','fontsize',fontsize)
    set(h_baselineplot,'PaperPositionMode','auto')
    print(h_baselineplot,'-depsc',[save_path filesep 'RMSE_baseline_plot.eps'])
    print(h_baselineplot,'-dpng',[save_path filesep 'RMSE_baseline_plot.png'])
end

% make a plot with what data is available where
no_correction_matrix = zeros([n_techniques+2 n_ifgs]);
for technique_counter=1:n_techniques+2
    no_correction_matrix(technique_counter,ix_no_correction{technique_counter}) = 1;
end
no_correction_matrix(end,((sum(no_correction_matrix(1:end-1,:),1)>=1)-no_correction_matrix(end,:))<0)=2;
no_correction_matrix(end,(no_correction_matrix(end,:)==1))=0;

% color based on technique availablilty
no_correction_matrix_color = no_correction_matrix;
for technique_counter=1:n_techniques+1
    no_correction_matrix_color(technique_counter,no_correction_matrix_color(technique_counter,:)==1)=technique_counter;
end
no_correction_matrix_color(end,no_correction_matrix_color(end,:)==2)=n_techniques+2;


% plotting the results
h_technique_matrix = figure('name','Technique numbers','position',[  200         267        1501         579]);
for temp_counter=1:length(legend_str)
   plot(1,1,'sq','color',current_colors(temp_counter,:),'markerfacecolor',current_colors(temp_counter,:))
   hold on 
end
plot(1,1,'sq','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5])
hold on
imagesc(1:n_ifgs,1:n_techniques+2,no_correction_matrix_color)
hold on
for technique_counter=1:n_techniques+3
   hold on
   plot([1-0.5 n_ifgs+0.5],[technique_counter technique_counter]-0.5,'k-')
end
for ifg_counter=1:n_ifgs+1
   hold on
   plot([ifg_counter ifg_counter]-0.5,[1-0.5 n_techniques+2+0.5],'k-')
end
axis equal
axis tight
axis ij
box on
title(['No correction estimated for marked squares'],'fontsize',fontsize)
xlabel('Interferogram number','fontsize',fontsize)
legend_str_temp =legend_str;
legend_str_temp{end+1} = ['Did not meet ' num2str(perc_threshold) '% threshold'];
legend('location','northoutside',legend_str_temp)
set(gca,'fontsize',fontsize,'ytick',[])
caxis([0 n_techniques+2])
colormap([1 1 1 ;current_colors ; 0.5 0.5 0.5])
set(h_technique_matrix,'PaperPositionMode','auto')
print(h_technique_matrix,'-depsc',[save_path filesep 'technique_matrix.eps'])
print(h_technique_matrix,'-dpng',[save_path filesep 'technique_matrix.png'])


% plotting the interferograms used in the RMSE computation for the
% different techniques
str = '';
while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
    str = input(['Plot the infividual interferogram comparison? [y/n] \n'],'s');
end
% only do this when requested
if strcmpi(str,'y')
    % getting the colarbar extends
    temp = reshape(data_option1(:,ix_ifgs,:),[],n_techniques+1);
    temp(isnan(temp(:,1)),:)=[];
    temp = sort(temp);
    ix_025perc = floor(size(temp,1)*0.025);
    if ix_025perc<=0
        ix_025perc=1;
    end
    ix_975perc = ceil(size(temp,1)*0.975);
    if ix_975perc<=0
        ix_975perc=size(temp,1);
    end
    cbar = [min([temp(ix_025perc,:)]) max([temp(ix_975perc,:)])];
    clear ix_025perc ix_975perc temp

    for technique_counter=1:n_techniques+1
        ps_plot(data_option1(:,ix_ifgs,technique_counter),0,cbar)   
        h_technique_plot = gcf;
        set(h_technique_plot,'name',save_str{technique_counter})

        % the first plot is the reference
        if technique_counter==1
            savename = [save_path filesep 'data_ref_' save_str{technique_counter}];
        % all these plots are techniques
        else
            savename = [save_path filesep 'data_' save_str{technique_counter}];
        end
        % saving the plots
        set(h_technique_plot,'PaperPositionMode','auto')
        print(h_technique_plot,'-depsc',[savename '.eps'])
        print(h_technique_plot,'-dpng',[savename '.png'])
    end
end


% saving the data used to make this plot
if exist([save_path filesep 'RMSE_output.mat'],'file')==2
    save([save_path filesep 'RMSE_output.mat'],'-append','scatterplot_mean','scatterplot_ref','scatterplot_ref_mean','units_flag','scatterplot','legend_str','RMSE_perc','ref_technique','other_techniques','perc_threshold','deramp_flag','reference_flag','crop_out_flag','crop_file_path','ix_no_correction','no_correction_matrix_color')
else
    save([save_path filesep 'RMSE_output.mat'],'scatterplot_mean','scatterplot','scatterplot_ref','scatterplot_ref_mean','units_flag','legend_str','RMSE_perc','ref_technique','other_techniques','perc_threshold','deramp_flag','reference_flag','crop_out_flag','crop_file_path','ix_no_correction','no_correction_matrix_color')
end


