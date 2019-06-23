function [] = aps_weather_model_InSAR(model_type)
% [] = aps_weather_model_InSAR(model_type)
% Goes to the InSAR data path and interpolates the weathermodel data to the InSAR
% grid. The tropospheric correction results are stored in the "tca2.mat" or 
% "tca_sb2.mat" file as the ph_tropo_era for ERA, and ph_tropo_merra for MERRA.
% The sign convention is defined such ph_after_corection = ph - ph_tropo_* is the phase corrected 
% the tropospheric signal. 
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
% By David Bekaert -- University of Leeds
%
% Modifications
% 10/2013   DB      Include conversion from zenith to slant delay 
% 10/2013   DB      Include processing structure different from stamps
% 11/2013   DB      Change the name to Hydrostratic as dry is incorrect
%                   naming
% 01/2014   DB      Fix in the saving of the tca filename.
% 05/2014   DB      Adding compatibility with non-stamps processing
%                   structures
% 06/2014   DB      Add check to make sure both master and slave have data.
% 02/2015   DB      Fix for the look angle in case its a single value
% 04/2016   DB      Convert code to be weather model independent and use
%                   modular approach
% 04/2016   DB      Forgot to remove call to old script
% 05/2016   DB      Include merra2
% 10/2016   DB 	    Include a fix inconsistent ifg matrix for stamps users who dropped ifgs
% 10/2016   DB      Change to aps_save command with append functionality.
% 11/2017   DB      Adding GACOS support, update the wording to inc angle
%                   as that the correct name

%% ERA-I, ERA5, WRF, MERRA1-2
% Filename suffix of the output files
wetoutfile = '_ZWD.xyz';   
hydroutfile = '_ZHD.xyz'; 
%% GACOS
% Filename suffix of the output files
outfile = '.ztd';

% error calling if needed
if nargin<1
    error('Give at least the model_type: era or merra')
end
% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed',1);
ll_matfile = getparm_aps('ll_matfile',1);
ifgday_matfile = getparm_aps('ifgday_matfile');

model_type = lower(model_type);
if strcmpi(model_type,'era')
    weather_model_datapath = getparm_aps('era_datapath',1);
elseif strcmpi(model_type,'merra') || strcmpi(model_type,'merra2')
    weather_model_datapath = getparm_aps('merra_datapath',1);
elseif strcmpi(model_type,'gacos')
    weather_model_datapath = getparm_aps('gacos_datapath',1);
else
    error('Not a supported model: either: ERA, MERRA, MERRA2, GACOS')
end 
lambda = getparm_aps('lambda',1)*100;                       % radar wavelength in cm
datestructure = 'yyyymmdd';                               % assumed date structure for era
inc_angle =  getparm_aps('look_angle',1);
% check if wavelength is specified
if isempty(lambda)
    error('Specify the wavelength, lambda is empty')
end

% loading the data
if strcmp(stamps_processed,'y')
   fprintf('Stamps processed structure \n')
   ps = load(ll_matfile);
   load psver
   dates = ps.day;
   lonlat = ps.lonlat;
   if ischar(inc_angle)==1
       inc_angle = load(inc_angle);
       inc_angle = inc_angle.la;
   end
   
   % getting the dropped ifgs
   drop_ifg_index = getparm('drop_ifg_index');
    % getting the parms file list from stamps to see the final ifg list
    if strcmp(getparm('small_baseline_flag'),'y')
        sb_flag = 1;
    else
        sb_flag = 0;
    end
    
    n_ifg = ps.n_ifg;
    % constructing the matrix with master and slave dates
    if sb_flag ==1
        % for SB
        ifg_number = [1:n_ifg]';
        ifgday_ix = ps.ifgday_ix;
        % removing those dropped interferograms
        ifgday_ix(drop_ifg_index,:) =[];
        ifg_number(drop_ifg_index)=[];

        % defining ix interferograms for which the delay needs to be computed
        ifgs_ix = [ifgday_ix ifg_number];
    else
        % slightly different for PS.
        date_slave_ix = [1:n_ifg]';
        ifg_number = [1:n_ifg]';

        % removing those interferograms that have been dropped
        date_slave_ix(drop_ifg_index)=[];
        ifg_number(drop_ifg_index)=[];

        % the master dates
        date_master_ix = repmat(ps.master_ix,size(date_slave_ix,1),1);

        % ix interferograms
        ifgs_ix = [date_master_ix date_slave_ix ifg_number];
    end
else
    psver = 2;
	
    % small baseline flag
    sb_flag = 0;
	    
    % loading lon lat information
    lonlat = load(ll_matfile);
    lonlat = lonlat.lonlat;
    
    % loading look angle information
    if ischar(inc_angle)==1
       inc_angle = load(inc_angle);
       inc_angle = inc_angle.la;
    end
    
    % getting the dates in jullian format
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
    dates = sort(dates);        % dates increasing with time
    
    % getting the ix position for the master and slave dates with respect
    % to the times
    date_master = datenum(num2str(ifgs_dates(:,1)),'yyyymmdd');
    date_slave = datenum(num2str(ifgs_dates(:,2)),'yyyymmdd');
    ifg_number = [1:size(date_master,1)]';
    
    drop_ifg_index=[];
    
    % dropping specific interferograms
    date_slave(drop_ifg_index)=[];
    date_master(drop_ifg_index)=[];
    ifg_number(drop_ifg_index)=[];
    for k=1:size(date_master,1)
        [date_master_ix(k,1)] = find(date_master(k,1)==dates);
        [date_slave_ix(k,1)] = find(date_slave(k,1)==dates);
    end
    
    % ix interferograms
    ifgs_ix = [date_master_ix date_slave_ix ifg_number];
    
    % number of total interferograms including dropped
    n_ifg = size(ifgs_dates,1);
    
    % defining ps.day information
    ps.day = dates;
end
n_dates = length(dates);
InSAR_datapath=['.'];
apsname = [InSAR_datapath filesep 'tca' num2str(psver) '.mat'];
apssbname = [InSAR_datapath filesep 'tca_sb' num2str(psver) '.mat'];



%% loading the weather model data
% initialisation 
d_wet = NaN([size(lonlat,1) n_dates]);          % these are the SAR estimated tropospheric delays for all data
d_hydro = NaN([size(lonlat,1) n_dates]);        % these are the SAR estimated tropospheric delays for all data
d_total = NaN([size(lonlat,1) n_dates]);        % these are the SAR estimated tropospheric delays for all data
flag_wet_hydro_used = 'n';

ix_no_weather_model_data = [];
counter = 0;

% looping over the dates 
for k=1:n_dates
    % getting the SAR data and convert it to a string
    date_str = datestr(ps.day(k,1),datestructure);

    % filenames
    model_filename_wet = [weather_model_datapath filesep date_str filesep date_str wetoutfile];
    model_filename_hydro = [weather_model_datapath filesep date_str filesep date_str hydroutfile];
    model_filename = [weather_model_datapath filesep date_str filesep date_str outfile];

    % checking if there is actual data for this date, if not just
    % leave NaN's in the matrix.
    
    if exist(model_filename_wet,'file') ==2
        flag_wet_hydro_used = 'y';
        % computing the dry delay
        [xyz_input,xyz_output] = load_weather_model_SAR(model_filename_hydro,lonlat);
        % saving the data which have not been interpolated
        d_hydro(:,k) = xyz_output(:,3);
        clear xyz_input xyz_output
        
        % computing the wet delays
        [xyz_input,xyz_output] = load_weather_model_SAR(model_filename_wet,lonlat);
         
        % saving the output data
        d_wet(:,k) = xyz_output(:,3);
        clear xyz_output
        counter = counter+1;
    elseif exist(model_filename,'file') ==2
        flag_wet_hydro_used = 'n';
        % this is the GACOS model file, will need to pass the model-type as
        % its grid-note 
        [xyz_input,xyz_output] = load_weather_model_SAR(model_filename,lonlat,[],model_type);
        % saving the output data
        d_total(:,k) = xyz_output(:,3);
        clear xyz_output
        counter = counter+1;

    else
        % rejected list of weather model images
       ix_no_weather_model_data = [ix_no_weather_model_data k]; 
    end
    clear model_filename_hydro model_filename_wet date_str
end
fprintf([num2str(counter) ' out of ' num2str(n_dates) ' SAR images have a tropospheric delay estimated \n'])




%% Computing the type of delay
if strcmpi(flag_wet_hydro_used,'y')
    d_total = d_hydro+d_wet;
end

%% Converting the Zenith delays to a slant delay
if size(inc_angle,2)>1 && size(inc_angle,1)==1
    inc_angle=inc_angle';
end
if size(inc_angle,2)==1
    inc_angle = repmat(inc_angle,1,size(d_total,2));
    if size(inc_angle,1)==1
        inc_angle = repmat(inc_angle,size(d_total,1),1);
    end
end

if strcmpi(flag_wet_hydro_used,'y')    
    d_hydro = d_hydro./cos(inc_angle);
    d_wet = d_wet./cos(inc_angle);
end
d_total = d_total./cos(inc_angle);


%% Converting the range delay to a phase delay
% converting to phase delay. 
% The sign convention is such that ph_corrected = ph_original - ph_tropo*
eval(['ph_SAR_' model_type '= -4*pi./lambda.*d_total;']);           % ph_SAR_era = -4*pi./lambda.*d_total;
if strcmpi(flag_wet_hydro_used,'y')
    eval(['ph_SAR_' model_type '_hydro= -4*pi./lambda.*d_hydro;']);     % ph_SAR_era_dry= -4*pi./lambda.*d_hydro;
    eval(['ph_SAR_' model_type '_wet= -4*pi./lambda.*d_wet;']);         % ph_SAR_era_wet= -4*pi./lambda.*d_wet;
end
clear d_total d_hydro d_wet

%% Computing the interferometric tropopsheric delays
% removing the dates for which there is no data.
if isempty(ix_no_weather_model_data)~=1
    for k=1:length(ix_no_weather_model_data)    
       % reject based on slave dates
       ix_ifg_reject = find(ix_no_weather_model_data(k)==ifgs_ix(:,2));
       ifgs_ix(ix_ifg_reject,:)=[];
       % reject based on master dates
       ix_ifg_reject = find(ix_no_weather_model_data(k)==ifgs_ix(:,1));
       ifgs_ix(ix_ifg_reject,:)=[];      
       clear ix_ifg_reject
    end
end

% computing the interferometric delay for each remaining interferogram
if isempty(ifgs_ix)
    fprintf(['Not enough ' upper(model_type) ' data to compute interferometric delays...\n'])
end
% initialize the ERA phase matrix for all interferograms, including those without correction.
eval(['ph_tropo_' model_type '= zeros([size(lonlat,1) n_ifg]);']);          % ph_tropo_era = zeros([size(lonlat,1) n_ifg]);
if strcmpi(flag_wet_hydro_used,'y')
    eval(['ph_tropo_' model_type '_hydro= zeros([size(lonlat,1) n_ifg]);']);    % ph_tropo_era_hydro = zeros([size(lonlat,1) n_ifg]);
    eval(['ph_tropo_' model_type '_wet= zeros([size(lonlat,1) n_ifg]);']);      % ph_tropo_era_wet = zeros([size(lonlat,1) n_ifg]);
end


n_ifg_kept = size(ifgs_ix,1);
for k=1:n_ifg_kept
    % add extra flag that requires both master and slave SAR delay to be
    % present otherwize, its still a sar delay!
    good_ifg_data = 1;
    % evaluate the following statement: sum(ph_SAR_era(:,ifgs_ix(k,1)))==0 |  sum(ph_SAR_era(:,ifgs_ix(k,2)))==0
    if eval(['sum(ph_SAR_' model_type '(:,ifgs_ix(' num2str(k) ',1)))==0 |  sum(ph_SAR_' model_type '(:,ifgs_ix(' num2str(k) ',2)))==0']);
        printf('Misses one of the SAR images for the interferometric delay \n')
        good_ifg_data = 0;
    end
    if good_ifg_data==1
        % ph_tropo_era(:,ifgs_ix(k,3)) = ph_SAR_era(:,ifgs_ix(k,1))-ph_SAR_era(:,ifgs_ix(k,2));
        eval(['ph_tropo_' model_type '(:,ifgs_ix(' num2str(k) ',3)) = ph_SAR_' model_type '(:,ifgs_ix(' num2str(k) ',1))-ph_SAR_' model_type '(:,ifgs_ix(' num2str(k) ',2));']);
        if strcmpi(flag_wet_hydro_used,'y')
            eval(['ph_tropo_' model_type '_hydro(:,ifgs_ix(' num2str(k) ',3)) = ph_SAR_' model_type '_hydro(:,ifgs_ix(' num2str(k) ',1))-ph_SAR_' model_type '_hydro(:,ifgs_ix(' num2str(k) ',2));']);
            eval(['ph_tropo_' model_type '_wet(:,ifgs_ix(' num2str(k) ',3)) = ph_SAR_' model_type '_wet(:,ifgs_ix(' num2str(k) ',1))-ph_SAR_' model_type '_wet(:,ifgs_ix(' num2str(k) ',2));']);
        end
    else
        % ph_tropo_era(:,ifgs_ix(k,3)) = NaN;
        eval(['ph_tropo_' model_type '(:,ifgs_ix(' num2str(k) ',3)) =NaN;']);
        if strcmpi(flag_wet_hydro_used,'y')
            eval(['ph_tropo_' model_type '_hydro(:,ifgs_ix(' num2str(k) ',3)) =NaN;']);
            eval(['ph_tropo_' model_type '_wet(:,ifgs_ix(' num2str(k) ',3)) =NaN;']);
        end
    end
end



if sb_flag==1
    save_name = apssbname;
else
    save_name = apsname;
end

if strcmpi(flag_wet_hydro_used,'y')
    eval(['aps_save(''' save_name ''',ph_tropo_' model_type ',ph_tropo_' model_type '_wet,ph_tropo_' model_type '_hydro);'])
else
    eval(['aps_save(''' save_name ''',ph_tropo_' model_type ');'])
end
