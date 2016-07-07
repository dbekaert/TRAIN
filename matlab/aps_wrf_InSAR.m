function [] = aps_wrf_InSAR()
% [] = aps_wrf_InSAR()
% Goes to the InSAR data path and interpolate the WRF-I result to the InSAR
% grid. The tropospheric correction results are stored in the "tca2.mat" or 
% "tca_sb2.mat" file as the ph_tropo_wrf variable.
% The sign convention is defined such ph_after_corection = ph - ph_tropo_wrf is the phase corrected 
% the tropospheric signal. 
%
% INPUTS:
%
%
% OPTIONAL INPUTS
% lambda                Radar wavelength [cm], by default c-band (5.6)
% datestructure         WRF date structure, by default this is assumed 'yyyymmdd'.
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
% 04/2015   DB      Bug correction for non-stamps processed case. ps
%                   variable definition of dates, definition of lonlat, look angle.

% Filename suffix of the output files
wetoutfile = '_ZWD.xyz';
hydroutfile = '_ZHD.xyz'; 

% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed');
ll_matfile = getparm_aps('ll_matfile');
wrf_datapath = getparm_aps('wrf_datapath');
lambda = getparm_aps('lambda')*100;                       % radar wavelength in cm
datestructure = 'yyyymmdd';                               % assumed date structure for wrf
% getting the dropped ifgs
drop_ifg_index = getparm('drop_ifg_index');

% loading the data
if strcmp(stamps_processed,'y')
   fprintf('Stamps processed structure \n')
   ps = load(ll_matfile);
   load psver
   dates = ps.day;
   lonlat = ps.lonlat;
   look_angle =  getparm_aps('look_angle');
   if ischar(look_angle)==1
       look_angle = load(look_angle);
       look_angle = look_angle.la;
   end
   
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
    
    lonlat = load(ll_matfile);
    lonlat = lonlat.lonlat;
    look_angle =  getparm_aps('look_angle');
    if ischar(look_angle)==1
       look_angle = load(look_angle);
       look_angle = look_angle.la;
    end
    
    % getting the dates in jullian format
    ifgday_matfile = getparm_aps('ifgday_matfile');
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
InSAR_datapath=['.' filesep ];
apsname = [InSAR_datapath filesep 'tca' num2str(psver) '.mat'];
apssbname = [InSAR_datapath filesep 'tca_sb' num2str(psver) '.mat'];




%% loading the WRF-I data
% initialisation 
d_wrf_wet = NaN([size(lonlat,1) n_dates]);          % these are the SAR estimated tropospheric delays for all data
d_wrf_dry = NaN([size(lonlat,1) n_dates]);          % these are the SAR estimated tropospheric delays for all data


ix_no_wrf = [];
counter = 0;
for k=1:n_dates
    % getting the SAR data and convert it to a string
    date_str = datestr(ps.day(k,1),datestructure);

    % WRF-I filename
    wrf_filename_wet = [wrf_datapath filesep date_str filesep date_str wetoutfile];
    wrf_filename_hydro = [wrf_datapath filesep date_str filesep date_str hydroutfile];

    % checking if there is actual data for this date, if not just
    % leave NaN's in the matrix.
    if exist(wrf_filename_wet,'file') ==2
        
        % computing the dry delay
        [xyz_input,xyz_output] = load_wrf_SAR(wrf_filename_hydro,lonlat);
        % saving the data which have not been interpolated
        
        d_wrf_dry(:,k) = xyz_output(:,3);
        clear xyz_input xyz_output
        
        % computing the wet delays
        [xyz_input,xyz_output] = load_wrf_SAR(wrf_filename_wet,lonlat);
         
        % saving the output data
        d_wrf_wet(:,k) = xyz_output(:,3);
        clear xyz_output
        counter = counter+1;

    else
        % rejected list of WRF-I images
       ix_no_wrf = [ix_no_wrf k]; 
    end
    clear wrf_filename_hydro wrf_filename_wet date_str
end
fprintf([num2str(counter) ' out of ' num2str(n_dates) ' SAR images have a tropospheric delay estimated \n'])


%% Computing the type of delay
d_wrf = d_wrf_dry+d_wrf_wet;

%% Converting the Zenith delays to a slant delay
if size(look_angle,2)>1 && size(look_angle,1)==1
    look_angle=look_angle';
end
if size(look_angle,2)==1
    look_angle = repmat(look_angle,1,size(d_wrf,2));
    if size(look_angle,1)==1
        look_angle = repmat(look_angle,size(d_wrf,1),1);
    end
end
    
d_wrf = d_wrf./cos(look_angle);
d_wrf_dry = d_wrf_dry./cos(look_angle);
d_wrf_wet = d_wrf_wet./cos(look_angle);


%% Converting the WRF-I delays to a phase delay
% converting to phase delay. 
% The sign convention is such that ph_corrected = ph_original - ph_wrf
ph_SAR_wrf = -4*pi./lambda.*d_wrf;
ph_SAR_wrf_dry= -4*pi./lambda.*d_wrf_dry;
ph_SAR_wrf_wet= -4*pi./lambda.*d_wrf_wet;


%% Computing the interferometric tropopsheric delays


% removing the dates for which there is no data.
if isempty(ix_no_wrf)~=1
    for k=1:length(ix_no_wrf)    
       % reject based on slave dates
       ix_ifg_reject = find(ix_no_wrf(k)==ifgs_ix(:,2));
       ifgs_ix(ix_ifg_reject,:)=[];
       % reject based on master dates
       ix_ifg_reject = find(ix_no_wrf(k)==ifgs_ix(:,1));
       ifgs_ix(ix_ifg_reject,:)=[];      
       clear ix_ifg_reject
    end
end

% computing the interferometric delay for each remaining interferogram
n_ifg_wrf = size(ifgs_ix,1);
if isempty(ifgs_ix)
    fprintf('Not enough WRF data to compute interferometric delays...\n')
end
% initialize the WRF phase matrix for all interferograms, including those without correction.
ph_tropo_wrf = zeros([size(lonlat,1) n_ifg]);
ph_tropo_wrf_hydro = zeros([size(lonlat,1) n_ifg]);
ph_tropo_wrf_wet = zeros([size(lonlat,1) n_ifg]);

for k=1:n_ifg_wrf
    ph_tropo_wrf(:,ifgs_ix(k,3)) = ph_SAR_wrf(:,ifgs_ix(k,1))-ph_SAR_wrf(:,ifgs_ix(k,2));
    ph_tropo_wrf_hydro(:,ifgs_ix(k,3)) = ph_SAR_wrf_dry(:,ifgs_ix(k,1))-ph_SAR_wrf_dry(:,ifgs_ix(k,2));
    ph_tropo_wrf_wet(:,ifgs_ix(k,3)) = ph_SAR_wrf_wet(:,ifgs_ix(k,1))-ph_SAR_wrf_wet(:,ifgs_ix(k,2));
end

if sb_flag==1
    if exist(apssbname,'file')==2
        save(apssbname,'-append','ph_tropo_wrf','ph_tropo_wrf_wet','ph_tropo_wrf_hydro')
    else
        save(apssbname,'ph_tropo_wrf','ph_tropo_wrf_wet','ph_tropo_wrf_hydro')       
    end
else
    if exist(apsname,'file')==2
        save(apsname,'-append','ph_tropo_wrf','ph_tropo_wrf_wet','ph_tropo_wrf_hydro')
    else
        save(apsname,'ph_tropo_wrf','ph_tropo_wrf_wet','ph_tropo_wrf_hydro')       
    end
end
