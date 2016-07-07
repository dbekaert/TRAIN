function [] = aps_modis_InSAR()
% [] = aps_modis_InSAR()
% Goes to the InSAR data path and interpolate the modis result to the InSAR
% grid. Skip those images that do not have modis data to correct for.
% Also only keep those modis tracks that cover a percentage of the InSAR
% track. The tropospheric correction results are stored in the "tca2.mat" or "tca_sb2.mat" file as the ph_tropo_modis variable.
% The sign convention is defined such ph_after_corection = ph - ph_tropo_modis is the phase corrected 
% the tropospheric signal. 
%
% OPTIONAL INPUTS
% meris_perc_coverage 	Minimum percentage of PS locations that should have no cloud coverage.
%                       By default this value is set to 80%.
% lambda                Radar wavelength [cm], by default c-band (5.6)
% datestructure         Modis date structure, by default this is assumed 'yyyymmdd'.
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
% By David Bekaert - May 2014
% University of Leeds
%
% modifications
% DB    08/2014     Output the parameters that are being lated from parms aps
% DB    08/2014     Include support for calibrated MODIS data
% DB    02/2016     Fix for the incidence angle in case a singel value is given



% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed',1);
ll_matfile = getparm_aps('ll_matfile',1);
modis_datapath = getparm_aps('modis_datapath',1);
modis_perc_coverage = getparm_aps('meris_perc_coverage',1); % percentage of PS locations that need to 
%                                                           have meris coverage. In case this is less
%                                                           the meris data is rejected.
lambda = getparm_aps('lambda',1)*100;                       % radar wavelength in cm
look_angle =  getparm_aps('look_angle',1);
ifgday_matfile = getparm_aps('ifgday_matfile',1);
modis_recalibrated = getparm_aps('modis_recalibrated',1);

% checking if its recalibrated modis data or not
if strcmpi(modis_recalibrated,'y')
    recal_str = 'recal_';
else
   recal_str = '';
end
% The modis files used for the estimation. 
% The no interp file is used to check the cloud coverage.
% the interpolated file is used as correction result.
modis_file_suffix_nointerp = [recal_str 'ZWD_nointerp.xyz'];
modis_file_suffix_interp = [recal_str 'ZWD_surf.xyz'];

datestructure = 'yyyymmdd';                               % assumed date structure for meris
% getting the dropped ifgs


% loading the data
if strcmp(stamps_processed,'y')
   fprintf('Stamps processed structure \n')
   drop_ifg_index = getparm('drop_ifg_index');
   
     % loading look angle information
    if ischar(look_angle)==1
       look_angle = load(look_angle);
       look_angle = look_angle.la;
    end

    % longitude, latitude and time information
   ps = load(ll_matfile);
   load psver
   dates = ps.day;
   lonlat = ps.lonlat;
   
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
        ifgs_ix_no_interp = ifgs_ix;

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
        ifgs_ix_no_interp = ifgs_ix;
    end

else
    
    sb_flag = 0;
    psver = 2;
    % loading lon lat information
    lonlat = load(ll_matfile);
    lonlat = lonlat.lonlat;
    
    % loading look angle information
    if ischar(look_angle)==1
       look_angle = load(look_angle);
       look_angle = look_angle.la;
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
    
    % dropping specific interferograms
    for k=1:size(date_master,1)
        [date_master_ix(k,1)] = find(date_master(k,1)==dates);
        [date_slave_ix(k,1)] = find(date_slave(k,1)==dates);
    end
    
    % ix interferograms
    ifgs_ix = [date_master_ix date_slave_ix ifg_number];
    
    % number of total interferograms including dropped
    n_ifg = size(ifgs_dates,1);
    
    fprintf('Check the InSAR dates, this has not been tested \n')
end


InSAR_datapath=['.' filesep ];
apsname = [InSAR_datapath 'tca' num2str(psver) '.mat'];
apssbname = [InSAR_datapath 'tca_sb' num2str(psver) '.mat'];



n_dates = length(dates);
n_points = size(lonlat,1);

%% loading the meris data
% initialisation 
d_modis = NaN([n_points n_dates]);          % these are the SAR estimated tropospheric delays for all data
d_modis_no_interp = NaN([n_points n_dates]);          % these are the SAR estimated tropospheric delays for all data


SAR_modis_perc = NaN([n_dates 1]);	   % Vector with the actual percentge of cloud cover. NaN for those were there is no meris data.
counter = 0;
counter_threshold = 0;
ix_no_modis =[];
ix_no_modis_coverage = [];
ix_no_modis_no_interp_file = [];

for k=1:n_dates
    % getting the SAR data and convert it to a string
    modis_date_str = datestr(dates(k,1),datestructure);

    % meris filename
    modis_filename_interp = [modis_datapath filesep modis_date_str filesep modis_date_str '_' modis_file_suffix_interp];
    modis_filename_nointerp = [modis_datapath filesep modis_date_str filesep modis_date_str '_' modis_file_suffix_nointerp];

    % checking if there is actual modis data for this date, if not just
    % leave NaN's in the matrix.
    if exist(modis_filename_interp,'file') ==2
        % computing the percentage of pixels that do not have data coverage
        if exist(modis_filename_nointerp,'file') ==2 
            [xyz_input_nointerp,xyz_output_nointerp] = load_meris_SAR(modis_filename_nointerp,lonlat);

            % saving the data which have not been interpolated
            d_modis_no_interp(:,k) = xyz_output_nointerp(:,3);

            SAR_modis_perc(k,1) = (1-sum(isnan(xyz_output_nointerp(:,3)))./n_points)*100;
            if (1-sum(isnan(xyz_output_nointerp(:,3)))./n_points)*100 < modis_perc_coverage
                % reject the modis acquisition as to much of it is covered
                % by clouds.
                keep_modis = 0;
                counter_threshold = counter_threshold+1;
            else
                % good modis acquisition
                keep_modis = 1;
            end
            clear xyz_input_nointerp xyz_output_nointerp
        else
            % no way to test if its a good modis track, keep it.
            % best is to keep all the outputs from the modis computation
            keep_modis = 1;
        end

        if keep_modis==1
            % load the interpolated modis data
            [xyz_input_interp,xyz_output_interp] = load_meris_SAR(modis_filename_interp,lonlat);
            clear xyz_input_interp 

            % saving the output data
            d_modis(:,k) = xyz_output_interp(:,3);
            clear xyz_output_interp

            % counting the number of SAR dates with modis data that have
            % enough coverage
            counter = counter +1;
        else
            % rejected list of modis images
           ix_no_modis = [ix_no_modis k]; 
        end
    else
        % rejected list of modis images
       ix_no_modis_no_interp_file = [ix_no_modis_no_interp_file k];
       ix_no_modis = [ix_no_modis k]; 
       ix_no_modis_coverage = [ix_no_modis_coverage k];
    end
    clear modis_filename modis_date_str

end

fprintf(['\n-----------------------------------------------------\n' num2str(counter) ' out of ' num2str(n_dates) ' SAR images have a tropospheric delay estimated \n'])
fprintf([num2str(counter_threshold) ' images did not meet the ' num2str(modis_perc_coverage) ' procent threshold. \n'])
% outputting a summary of the modis data:
no_modis_data = repmat(' ',n_dates,13);
no_modis_data(ix_no_modis_coverage,:)=repmat('no modis data',length(ix_no_modis_coverage),1);
output_str = [datestr(dates,datestructure) repmat(':  ',n_dates,1) num2str(round(SAR_modis_perc*10)/10) repmat(' procent     ',n_dates,1) no_modis_data repmat('\n',n_dates,1)];
for k=1:n_dates
    fprintf([output_str(k,:)])
end 
fprintf(['\n-----------------------------------------------------\n'])


%% Converting the Zenith delays to a slant delay
if size(look_angle,2)>1 && size(look_angle,1)==1
    look_angle=look_angle';
end
if size(look_angle,2)==1 && size(look_angle,2)~=1       % fix for single angles
    look_angle = repmat(look_angle,1,size(d_modis,2));
end
d_modis = d_modis./cos(look_angle);
d_modis_no_interp = d_modis_no_interp./cos(look_angle);


%% Converting the modis delays to a phase delay
% converting to phase delay. 
% The sign convention is such that ph_corrected = ph_original - ph_modis
ph_SAR_modis = -4*pi./lambda.*d_modis;
ph_SAR_modis_no_interp = -4*pi./lambda.*d_modis_no_interp;
clear d_modis d_modis_no_interp


%% Computing the interferometric tropopsheric delays
% removing the dates for which there is no data.
if isempty(ix_no_modis)~=1
    for k=1:length(ix_no_modis)    
       % reject based on slave dates
       ix_ifg_reject = find(ix_no_modis(k)==ifgs_ix(:,2));
       ifgs_ix(ix_ifg_reject,:)=[];
       % reject based on master dates
       ix_ifg_reject = find(ix_no_modis(k)==ifgs_ix(:,1));
       ifgs_ix(ix_ifg_reject,:)=[];      
       clear ix_ifg_reject
    end
end

% removing the dates for which there is no data.
if isempty(ix_no_modis_no_interp_file)~=1
    for k=1:length(ix_no_modis_no_interp_file)    
       % reject based on slave dates
       ix_ifg_reject_no_interp = find(ix_no_modis_no_interp_file(k)==ifgs_ix_no_interp(:,2));
       ifgs_ix_no_interp(ix_ifg_reject_no_interp,:)=[];
       % reject based on master dates
       ix_ifg_reject_no_interp = find(ix_no_modis_no_interp_file(k)==ifgs_ix_no_interp(:,1));
       ifgs_ix_no_interp(ix_ifg_reject_no_interp,:)=[];      
       clear ix_ifg_reject_no_interp
    end
end

% computing the interferometric delay for each remaining interferogram
n_ifg_modis = size(ifgs_ix,1);
if isempty(ifgs_ix)
    fprintf('Not enough modis data to compute interferometric delays...\n')
end
% initialize the modis phase matrix for all interferograms, including those
% without correction.
ph_tropo_modis = zeros([n_points n_ifg]);
ph_tropo_modis_no_interp = NaN([n_points n_ifg]);
for k=1:n_ifg_modis
    ph_tropo_modis(:,ifgs_ix(k,3)) = ph_SAR_modis(:,ifgs_ix(k,1))-ph_SAR_modis(:,ifgs_ix(k,2));
end
n_ifg_modis_no_interp = size(ifgs_ix_no_interp,1);
for k=1:n_ifg_modis_no_interp
    ph_tropo_modis_no_interp(:,ifgs_ix_no_interp(k,3)) = ph_SAR_modis_no_interp(:,ifgs_ix_no_interp(k,1))-ph_SAR_modis_no_interp(:,ifgs_ix_no_interp(k,2));
end



% check if the data needs to be saved as recalibrated modis or not
if strcmpi(modis_recalibrated,'y')          % recalibrated modis
    ph_tropo_modis_recal = ph_tropo_modis;
    ph_tropo_modis_no_interp_recal = ph_tropo_modis_no_interp;
    clear ph_tropo_modis_no_interp ph_tropo_modis
    if sb_flag==1
        if exist(apssbname,'file')==2
            save(apssbname,'-append','ph_tropo_modis_recal','ph_tropo_modis_no_interp_recal')
        else
            save(apssbname,'ph_tropo_modis_recal','ph_tropo_modis_no_interp_recal')        
        end
    else
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis_recal','ph_tropo_modis_no_interp_recal')
        else
            save(apsname,'ph_tropo_modis_recal','ph_tropo_modis_no_interp_recal')        
        end
    end
else                                        % not prior calibrated
    if sb_flag==1
        if exist(apssbname,'file')==2
            save(apssbname,'-append','ph_tropo_modis','ph_tropo_modis_no_interp')
        else
            save(apssbname,'ph_tropo_modis','ph_tropo_modis_no_interp')        
        end
    else
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis','ph_tropo_modis_no_interp')
        else
            save(apsname,'ph_tropo_modis','ph_tropo_modis_no_interp')        
        end
    end
end

% output cloud information
modis.SAR_no_cloud_perc = SAR_modis_perc;
modis.ix_no_data = ix_no_modis_coverage;
if exist('tca_support.mat','file')==2
    save('tca_support.mat','-append','modis')
else
    save('tca_support.mat','modis')        
end
