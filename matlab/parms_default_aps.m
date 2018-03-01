function parms_default_aps()
% lwts_parms_default set parms to default value if not already set
%
% Based on script by Andy Hooper (StaMPS)
% Modified for the aps toolbox by David Bekaert - University of Leeds - 2013
% This script allows for non-StaMPS structured processed data.
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
% 04/2013   DB:     Convert column values to a row when logging the changes
% 10/2013   DB:     Include ERA-I options
% 10/2013   DB:     Fixed flag for StaMPS SB option
% 10/2013   DB:     Fix bug for UTC time retrieval of master.res file
% 11/2013   DB:     Including WRF in parameter list
% 02/2014   DB:     Introduce power_law_ridge_constraint to use ridge
%                   information to constrain the patch definition
% 03/2014   DB:     Change default ERA-I to ECMWF
% 03/2014   DB:     Add auto selection of merid grid for StamPS users
% 05/2014   DB:     Add MODIS support
% 05/2014   DB:     replace variable names to spectrometer as its the same for modis and meris
% 05/2014   DB:     Fix error in default look angle value
% 07/2014   DB:     Add the modis calibration factor
% 08/2014   DB:     Add an extra powerlaw_kept flag
% 08/2014   DB:     Include support for gamma processed data UTC retrieval
% 08/2014   DB:     Add option to use recalibrated MODIS data
% 01/2016   DB:     Include powerlaw plane flag
% 04/2016   DB:     Adding MERRA support
% 02/2018   KM:     Adding NARR support

parmfile='parms_aps.mat';
parent_flag=0;

if exist(['.' filesep parmfile],'file')
    parms=load(parmfile);
elseif exist(['..' filesep parmfile],'file')
    parmfile=['..' filesep parmfile];
    parms=load(parmfile);
    parent_flag=1;
else
    parms=struct('Created',date);
    parms.small_baseline_flag='n'; 
end


parmfields_before=fieldnames(parms);
num_fields=size(parmfields_before,1);

warning_message = 0;
%% StaMPS specific
if exist('psver.mat')~=2
   stamps_flag = 0; 
   parms.stamps_processed = 'n';
else
    load psver
    psname = ['ps',num2str(psver)];
    ps = load(psname);
    stamps_flag = 1; 
    parms.stamps_processed = 'y';
end
%% regarding sounding data
if ~isfield(parms,'sounding_data')
    % use sounding data when 1
    parms.sounding_data='n';
end
if ~isfield(parms,'sounding_dir')
    % the sounding data dir
    parms.sounding_dir=[pwd filesep 'sounding_data'];
end
if ~isfield(parms,'sounding_time_stamp')
    % The sounding aquisition times
    parms.sounding_time_stamp=['00' ; '12'];
end
if ~isfield(parms,'sounding_start_date')
    % The sounding aquisition times. When [] use the first date.
    parms.sounding_start_date=[];
end
if ~isfield(parms,'sounding_end_date')
    % The sounding aquisition times. When [] use last date.
    parms.sounding_end_date=[];
end
if ~isfield(parms,'sounding_h0')
    % The maximum height at which the delay is computed in km. When 0 estimate
    % this height from the sounding data.
    parms.sounding_h0=0;
end
if ~isfield(parms,'sounding_error_promp')
    % sounding errors being outputed but keep on running script by setting
    % to 0.
    parms.sounding_error_promp='n';
end
if ~isfield(parms,'sounding_h_alpha_thres')
    % lower height range which is used to compute the powerlaw decay
    % coefficient, given in km
    parms.sounding_h_alpha_thres=4;
end

if ~isfield(parms,'sounding_months')
    % monthly interval for which the sensitivity is is evaluated
    parms.sounding_months=1;
end

if ~isfield(parms,'sounding_ifg_dates')
    % do sensitivity analysis only for SAR dates
    parms.sounding_ifg_dates='n';
end

if ~isfield(parms,'sounding_sensitivity')
    % perform sensitivity analysis
    parms.sounding_sensitivity='n';
end

%% satellite specific
if ~isfield(parms,'lambda')
    if stamps_flag==1
        % radar wavelength in [m]
        lambda_file = 'lambda.1.in';
        if exist(lambda_file,'file')~=2
            lambda_file = ['..' filesep lambda_file];
            if exist(lambda_file,'file')~=2
                lambda_file = ['..' filesep lambda_file];
            end
        end
        if exist(lambda_file,'file')==2
            parms.lambda = load(lambda_file);
        else
            fprintf('Set lambda manual in parms_aps! \n')
        end
    else
        parms.lambda = 0.0562;
        warning_message = 1;
    end
end
if ~isfield(parms,'look_angle')
    if stamps_flag==1
        % mean look angle in [deg]
        look_angle = [pwd filesep 'la'  num2str(psver) '.mat'];
        parms.look_angle = look_angle;
    else
        parms.look_angle = 21*pi/180;
        warning_message = 1;
    end
end
if ~isfield(parms,'heading')
    if stamps_flag==1
        % mean heading in [deg]
        heading_file = 'heading.1.in';
        if exist(heading_file,'file')~=2
            heading_file = ['..' filesep heading_file];
            if exist(heading_file,'file')~=2
                heading_file = ['..' filesep heading_file];
            end
        end
        if exist(heading_file,'file')==2
            parms.heading = load(heading_file);
        else
            fprintf('Set heading manual in parms_aps! \n')
        end
    else
        parms.heading = [];
        warning_message = 1;
    end
end
if ~isfield(parms,'UTC_sat')
    if stamps_flag==1
        % Satellite heading as a string in 'HH:MM'
        if strcmpi(getparm('insar_processor'),'gamma')
            load ps2.mat
            command_str = ['ls ' datestr(master_day,'yyyymmdd') '*' filesep datestr(master_day,'yyyymmdd') '*slc.par > temp.master.res'];
            [a,b] = system(command_str);
            fid=fopen('temp.master.res');
            
            if fid>0
                % getting only one master file back
                 master_files = textscan(fid,'%s');
                 master_file = master_files{1}(1);
                 
                 % dumping the time information for this file in a UTC file.
                 command_str = ['echo `grep center\_time\: ' char(master_file) '` > temp.utc'];
                 
                 [a,b] = system(command_str);
                 
                 if exist('temp.utc')~=2
                    warning_message = 1;
                    UTC_sat=[];
                 else
                     filetext = strtrim(fileread('temp.utc'));
                     ix = find(filetext==':');
                     [a,b] = system('rm temp.utc');
                     UTC_sat = str2num(filetext(ix(1)+1:end-1));
                     UTC_hr = floor(UTC_sat/60/60);
                     UTC_min = floor((UTC_sat-UTC_hr*60)/60);
                     UTC_hr = num2str(UTC_hr);
                     UTC_min = num2str(UTC_min);
                     if length(UTC_min)==1
                         UTC_min = ['0' UTC_min];
                     end
                     if length(UTC_hr)==1
                         UTC_hr = ['0' UTC_hr];
                     end
                     UTC_sat = [UTC_hr ':' UTC_min];
                 end
                 clear a b 
            end
            fclose(fid);
        else
            if exist('master.res')~=2
                if exist ('../master.res')==2
                    command_str = 'echo `grep First\_pixel\_azimuth\_time\ \(UTC\)\: ../master.res` > temp.utc';
                    [a,b] = system(command_str);
                end
            else
                command_str = 'echo `grep First\_pixel\_azimuth\_time\ \(UTC\)\: master.res` > temp.utc';
                [a,b] = system(command_str);
            end
            if exist('temp.utc')~=2
                warning_message = 1;
                UTC_sat=[];
            else
                filetext = fileread('temp.utc');
                ix = find(filetext==':');
                [a,b] = system('rm temp.utc');
                UTC_sat = filetext(ix(2)-2:ix(3)-1);
                clear a b 
            end
        end
    else
        warning_message = 1;
        UTC_sat=[];
    end
    parms.UTC_sat = UTC_sat;
end

%% regarding the aps estimation method
if ~isfield(parms,'non_defo_flag')
    parms.non_defo_flag='n';
end

%% region specific
if ~isfield(parms,'region_res')
    % dem_resolution in degrees
    parms.region_res = 0.008333;
end

if ~isfield(parms,'region_lon_range')
    % dem_resolution in degrees
    if stamps_flag==1
        temp = load([pwd filesep 'ps'  num2str(psver) '.mat']);
        lonlat_temp = temp.lonlat;
        clear temp
        if ~isempty(lonlat_temp)
            parms.region_lon_range = [floor(min(lonlat_temp(:,1)))-5*parms.region_res ceil(max(lonlat_temp(:,1)))+5*parms.region_res];
            clear lonlat_temp
        else
            parms.region_lon_range = [];
        end
    else
        parms.region_lon_range = [];
    end
end
if ~isfield(parms,'region_lat_range')
    % dem_resolution in degrees
    if stamps_flag==1
        temp = load([pwd filesep 'ps'  num2str(psver) '.mat']);
        lonlat_temp = temp.lonlat;
        clear temp
        if ~isempty(lonlat_temp)
            parms.region_lat_range = [floor(min(lonlat_temp(:,2)))-5*parms.region_res ceil(max(lonlat_temp(:,2)))+5*parms.region_res];
            clear lonlat_temp
        else
            parms.region_lat_range = [];
        end
    else
        parms.region_lat_range = [];
    end
end


if ~isfield(parms,'demfile')
    % dem_resolution in degrees
    parms.demfile = [pwd filesep 'dummy.dem'];
end
%% regarding spectrometer correction:
if ~isfield(parms,'spectrometer_scaleheight')
    % Scale height in m.
    parms.spectrometer_scaleheight = 8340;
end
if ~isfield(parms,'spectrometer_PIconversion')
    % PI conversion factor 
    parms.spectrometer_PIconversion = 6.2;
end

%% regarding the meris correction

if ~isfield(parms,'dem_null')
    % dem_null values
    parms.dem_null = -32768;
end


if ~isfield(parms,'meris_perc_coverage')
    % dem_resolution in degrees
    parms.meris_perc_coverage = 80;
end
if ~isfield(parms,'meris_datapath')
    % meris data path
    parms.meris_datapath = [];
end


%% regarding MODIS
if ~isfield(parms,'modis_datapath')
    % meris data path
    parms.modis_datapath = [];
end
% modis calibration factor from:
% Li, Z., Fielding, E. J., Cross, P., & Preusker, R. (2009). 
% Advanced InSAR atmospheric correction: MERIS/MODIS combination and 
% stacked water vapour models. International Journal of Remote Sensing, 30(13), 
% 3343-3363. doi: 10.1080/01431160802562172 
if ~isfield(parms,'modis_calibration')
    % meris data path
    parms.modis_calibration = 0.95;
end
if ~isfield(parms,'modis_recalibrated')
    % meris data path
    parms.modis_recalibrated = 'n';
end


%% regarding ERA-I
if ~isfield(parms,'era_datapath')
    % ERA-I data path
    parms.era_datapath = [];
end

if ~isfield(parms,'era_data_type')
    % ERA-I data website (BADC or ECMWF)
    parms.era_data_type ='ECMWF';
end

%% regarding gacos
if ~isfield(parms,'gacos_datapath')
    % gacos data path
    parms.gacos_datapath = [];
end


%% regarding MERRA
if ~isfield(parms,'merra_datapath')
    % MERRA data path
    parms.merra_datapath = [];
end

if ~isfield(parms,'era_data_type')
    % ERA-I data website (BADC or ECMWF)
    parms.era_data_type ='ECMWF';
end

%% regarding NARR
if ~isfield(parms,'narr_datapath')
    % NARR data path
    parms.narr_datapath = [];
end


%% regarding WRF
if ~isfield(parms,'wrf_datapath')
    % WRF data path
    parms.wrf_datapath = [];
end


%% regarding the powerlaw correction
if ~isfield(parms,'powerlaw_DEM_corr')
    % Powerlay decay coefficient.
    parms.powerlaw_DEM_corr = 'n';
end

if ~isfield(parms,'powerlaw_h0')
    % absolute height in km at which net-delay has reduced to zero. Default
    % this is 15m.
    parms.powerlaw_h0 = 10; 
end
if ~isfield(parms,'powerlaw_n_patches')
    parms.powerlaw_n_patches=50;
end
if ~isfield(parms,'powerlaw_alpha')
    % Powerlay decay coefficient.
    parms.powerlaw_alpha = 1.6;
end

if ~isfield(parms,'powerlaw_xy_res')
    % Grid resolution in m for the local grid.
    parms.powerlaw_xy_res = [30 30];
end
if ~isfield(parms,'powerlaw_patch_overlap')
    % percentage of the overlapping patches.
    parms.powerlaw_patch_overlap = 50;
end
if ~isfield(parms,'powerlaw_all_bands')
    % which band to use
    parms.powerlaw_all_bands = 'y';
end
if ~isfield(parms,'powerlaw_spatial_bands')
    % spatial bands filters.
    spatial_bands = [2000 4000
        4000    8000
        8000    16000
        16000   32000
        32000   64000
        64000  128000];
    parms.powerlaw_spatial_bands = spatial_bands;
end

if ~isfield(parms,'powerlaw_ridge_constraint')
    parms.powerlaw_ridge_constraint = 'n';
end

if ~isfield(parms,'powerlaw_kept')
    parms.powerlaw_kept = 0;
end

if ~isfield(parms,'powerlaw_plane_mode')
    parms.powerlaw_plane_mode = 'y';
end

%% ifg specific
if ~isfield(parms,'crop_flag')
    parms.crop_flag='n';
end


% path to themat file containing the interferogram dates
if ~isfield(parms,'ifgday_matfile')
    if stamps_flag==1
        parms.ifgday_matfile = [pwd filesep 'ps'  num2str(psver) '.mat'];
    else
        parms.ifgday_matfile=[];
        warning_message = 1;
    end

end


if ~isfield(parms,'save_folder_name')
    % save folder of the results
    parms.save_folder_name ='aps_estimation';
end

if ~isfield(parms,'drop_ifg_index')
    if stamps_flag==1
        % list of interferograms to be dropped out of the processing
        parms.drop_ifg_index = getparm('drop_ifg_index');
    else
        parms.drop_ifg_index = [];
    end
end

if ~isfield(parms,'phuw_matfile')
    if stamps_flag==1
        % list of interferograms to be dropped out of the processing
        small_baseline_flag = getparm('small_baseline_flag');
        parms.small_baseline_flag = small_baseline_flag;
        load psver
        if strcmp(small_baseline_flag,'y')
           % this is small baselines
           parms.phuw_matfile = [pwd filesep 'phuw_sb'  num2str(psver) '.mat'];
        else
           % this is single master approach
           parms.phuw_matfile = [pwd filesep 'phuw' num2str(psver) '.mat'];
        end
    else
       parms.phuw_matfile = [pwd filesep 'phuw.mat'];
       warning_message = 1;
    end
end
if ~isfield(parms,'hgt_matfile')
    if stamps_flag==1
        load psver
        % list of interferograms to be dropped out of the processing
        parms.hgt_matfile = [pwd filesep 'hgt'  num2str(psver) '.mat'];
    else
        parms.hgt_matfile = [pwd filesep 'hgt.mat'];
    end
end
if ~isfield(parms,'ll_matfile')
    if stamps_flag==1
        % list of interferograms to be dropped out of the processing
        parms.ll_matfile = [pwd filesep 'ps'  num2str(psver) '.mat'];
    else
        parms.ll_matfile = [pwd filesep 'll.mat'];
        warning_message = 1;
    end
end
if ~isfield(parms,'bperp_matfile')
    if stamps_flag==1
        % list of interferograms to be dropped out of the processing
        parms.bperp_matfile = [pwd filesep 'ps'  num2str(psver) '.mat'];
    else
        parms.bperp_matfile = [pwd filesep 'bperp.mat'];
        warning_message = 1;
    end
end


%%
parmfields=fieldnames(parms);
if size(parmfields,1)~=num_fields
    try
        save(parmfile,'-struct','parms')
        for i=1:size(parmfields,1)
            if isempty(strmatch(parmfields{i},parmfields_before))
               parmname=parmfields{i};
               value=getfield(parms,parmname);
               if isempty(value)
                   value='[]';
               end
               if isnumeric(value)               
                   if size(value,1)>1                   
                       value_str = '[';
                       for ll=1:size(value,1)
                           if ll==size(value,1)
                               value_str = [value_str num2str(value(ll,:)) ];
                           else
                                value_str = [value_str num2str(value(ll,:)) '; '];
                           end
                       end
                       value_str = [value_str ']'];
                   else
                       value_str = num2str(value);
                   end
                   logit([parmname,' = ',value_str],'aps.log',parent_flag);
                   clear value_str
               else
                   if size(value,1)>1
                      value_str = '[';
                      for ll=1:size(value,1)
                           if ll==size(value,1)
                               value_str = [value_str num2str(value(ll,:)) ];
                           else
                                value_str = [value_str num2str(value(ll,:)) '; '];
                           end
                       end
                       value_str = [value_str ']'];
                   else
                       value_str = value;
                   end
                   logit([parmname,' = ',value_str],'aps.log',parent_flag);
                   clear value_str
               end
            end
        end
        
    catch
        fprintf('Warning: missing parameters could not be updated (no write access)\n')
    end
end


if stamps_flag==0 && warning_message == 1
   fprintf(['\n ** Please check as default values were assumed for: \n lambda, heading, phuw_matfile, bperp_matfile, hgt_matfile, ll_matfile, demfile, region_lat_range, region_lon_range, spectrometer_scaleheight, spectrometer_PIconversion and dem_null \n\n\n']) 
end



 
