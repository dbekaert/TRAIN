function [] = aps_meris_InSAR()
% [] = aps_meris_InSAR()
% Goes to the InSAR data path and interpolate the meris result to the InSAR
% grid. Skip those images that do not have meris data to correct for.
% Also only keep those meris tracks that cover a percentage of the InSAR
% track. The tropospheric correction results are stored in the "tca2.mat" or "tca_sb2.mat" file as the ph_tropo_meris variable.
% The sign convention is defined such ph_after_corection = ph - ph_tropo_meris is the phase corrected 
% the tropospheric signal. Note that it only uses a wet delay. The dry
% delay is not ideal, as it is based on a exponential decay and scale
% height and it therefore not included. Use the ERA, WRF estiamted
% hydrostatic delays instead.
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
%
% By David Bekaert - September 2013
% University of Leeds
%
% modifications
% 10/2013   DB      Integrating other processing structure than stamps
% 05/2014   DB      Include variable name for non-stamps processing.
% 05/2015   DB      Output the information about percentage of cloud cover
% 06/2015   DB      flag the parameters that are used.

% The meris files used for the estimation. 
% The no interp file is used to check the cloud coverage.
% the interpolated file is used as correction result.
meris_file_suffix_nointerp = 'SWD_nointerp.xyz';
meris_file_suffix_interp = 'SWD_surf.xyz';


% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed',1);
ll_matfile = getparm_aps('ll_matfile',1);
meris_datapath = getparm_aps('meris_datapath',1);
meris_perc_coverage = getparm_aps('meris_perc_coverage',1); % percentage of PS locations that need to 
%                                                           have meris coverage. In case this is less
%                                                           the meris data is rejected.
lambda = getparm_aps('lambda',1)*100;                     % radar wavelength in cm
datestructure = 'yyyymmdd';                               % assumed date structure for meris
% getting the dropped ifgs


% loading the data
if strcmp(stamps_processed,'y')
   fprintf('Stamps processed structure \n')
   drop_ifg_index = getparm('drop_ifg_index');

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
    psver = 2;
    sb_flag = 0;
    
    % loading lon lat information
    lonlat = load(ll_matfile);
    lonlat = lonlat.lonlat;
    

    
    
   
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
d_meris = NaN([n_points n_dates]);          % these are the SAR estimated tropospheric delays for all data
d_meris_no_interp = NaN([n_points n_dates]);          % these are the SAR estimated tropospheric delays for all data


SAR_meris_perc = NaN([n_dates 1]);	   % Vector with the actual percentge of cloud cover. NaN for those were there is no meris data.
counter = 0;
counter_threshold = 0;
ix_no_meris =[];
ix_no_meris_coverage = [];
ix_no_meris_no_interp_file = [];
for k=1:n_dates
    % getting the SAR data and convert it to a string
    meris_date_str = datestr(dates(k,1),datestructure);

    % meris filename
    meris_filename_interp = [meris_datapath filesep meris_date_str filesep meris_date_str '_' meris_file_suffix_interp];
    meris_filename_nointerp = [meris_datapath filesep meris_date_str filesep meris_date_str '_' meris_file_suffix_nointerp];

    % checking if there is actual meris data for this date, if not just
    % leave NaN's in the matrix.
    if exist(meris_filename_interp,'file') ==2
        % computing the percentage of pixels that do not have data coverage
        if exist(meris_filename_nointerp,'file') ==2 
            [xyz_input_nointerp,xyz_output_nointerp] = load_meris_SAR(meris_filename_nointerp,lonlat);

            % saving the data which have not been interpolated
            d_meris_no_interp(:,k) = xyz_output_nointerp(:,3);

            SAR_meris_perc(k,1) = (1-sum(isnan(xyz_output_nointerp(:,3)))./n_points)*100;
            if (1-sum(isnan(xyz_output_nointerp(:,3)))./n_points)*100 < meris_perc_coverage
                % reject the meris acquisition as to much of it is covered
                % by clouds.
                keep_meris = 0;
                counter_threshold = counter_threshold+1;
            else
                % good meris acquisition
                keep_meris = 1;
            end
            clear xyz_input_nointerp xyz_output_nointerp
        else
            % no way to test if its a good meris track, keep it.
            % best is to keep all the outputs from the meris computation
            keep_meris = 1;
        end

        if keep_meris==1
            % load the interpolated meris data
            [xyz_input_interp,xyz_output_interp] = load_meris_SAR(meris_filename_interp,lonlat);
            clear xyz_input_interp 

            % saving the output data
            d_meris(:,k) = xyz_output_interp(:,3);
            clear xyz_output_interp

            % counting the number of SAR dates with meris data that have
            % enough coverage
            counter = counter +1;
        else
            % rejected list of meris images
           ix_no_meris = [ix_no_meris k]; 
        end
    else
        % rejected list of meris images
       ix_no_meris_no_interp_file = [ix_no_meris_no_interp_file k];
       ix_no_meris = [ix_no_meris k]; 
       ix_no_meris_coverage = [ix_no_meris_coverage k];
    end
    clear meris_filename meris_date_str

end

fprintf([num2str(counter) ' out of ' num2str(n_dates) ' SAR images have a tropospheric delay estimated \n'])
fprintf([num2str(counter_threshold) ' images did not meet the ' num2str(meris_perc_coverage) ' procent threshold. \n'])
% outputting a summary of the meris data:
no_meris_data = repmat(' ',n_dates,13);
no_meris_data(ix_no_meris_coverage,:)=repmat('no meris data',length(ix_no_meris_coverage),1);

output_str = [datestr(dates,datestructure) repmat(':  ',n_dates,1) num2str(round(SAR_meris_perc*10)/10) repmat(' procent     ',n_dates,1) no_meris_data repmat('\n  ',n_dates,1)];
for k=1:n_dates
    fprintf([output_str(k,:)])
end 
fprintf([num2str(n_dates-(counter+counter_threshold)) ' images did not have meris data \n'])



%% Converting the meris delays to a phase delay
% converting to phase delay. 
% The sign convention is such that ph_corrected = ph_original - ph_meris
ph_SAR_meris = -4*pi./lambda.*d_meris;
ph_SAR_meris_no_interp = -4*pi./lambda.*d_meris_no_interp;
clear d_meris d_meris_no_interp


%% Computing the interferometric tropopsheric delays
% removing the dates for which there is no data.
if isempty(ix_no_meris)~=1
    for k=1:length(ix_no_meris)    
       % reject based on slave dates
       ix_ifg_reject = find(ix_no_meris(k)==ifgs_ix(:,2));
       ifgs_ix(ix_ifg_reject,:)=[];
       % reject based on master dates
       ix_ifg_reject = find(ix_no_meris(k)==ifgs_ix(:,1));
       ifgs_ix(ix_ifg_reject,:)=[];      
       clear ix_ifg_reject
    end
end

% removing the dates for which there is no data.
if isempty(ix_no_meris_no_interp_file)~=1
    for k=1:length(ix_no_meris_no_interp_file)    
       % reject based on slave dates
       ix_ifg_reject_no_interp = find(ix_no_meris_no_interp_file(k)==ifgs_ix_no_interp(:,2));
       ifgs_ix_no_interp(ix_ifg_reject_no_interp,:)=[];
       % reject based on master dates
       ix_ifg_reject_no_interp = find(ix_no_meris_no_interp_file(k)==ifgs_ix_no_interp(:,1));
       ifgs_ix_no_interp(ix_ifg_reject_no_interp,:)=[];      
       clear ix_ifg_reject_no_interp
    end
end


% computing the interferometric delay for each remaining interferogram
n_ifg_meris = size(ifgs_ix,1);
if isempty(ifgs_ix)
    fprintf('Not enough meris data to compute interferometric delays...\n')
end
% initialize the meris phase matrix for all interferograms, including those
% without correction.
ph_tropo_meris = zeros([n_points n_ifg]);
ph_tropo_meris_no_interp = NaN([n_points n_ifg]);
for k=1:n_ifg_meris
    ph_tropo_meris(:,ifgs_ix(k,3)) = ph_SAR_meris(:,ifgs_ix(k,1))-ph_SAR_meris(:,ifgs_ix(k,2));
end
n_ifg_meris_no_interp = size(ifgs_ix_no_interp,1);
for k=1:n_ifg_meris_no_interp
    ph_tropo_meris_no_interp(:,ifgs_ix_no_interp(k,3)) = ph_SAR_meris_no_interp(:,ifgs_ix_no_interp(k,1))-ph_SAR_meris_no_interp(:,ifgs_ix_no_interp(k,2));
end

if sb_flag==1
    if exist(apssbname,'file')==2
        save(apssbname,'-append','ph_tropo_meris','ph_tropo_meris_no_interp')
    else
        save(apssbname,'ph_tropo_meris','ph_tropo_meris_no_interp')        
    end
else
    if exist(apsname,'file')==2
        save(apsname,'-append','ph_tropo_meris','ph_tropo_meris_no_interp')
    else
        save(apsname,'ph_tropo_meris','ph_tropo_meris_no_interp')        
    end
end


% output cloud information
meris.SAR_no_cloud_perc = SAR_meris_perc;
meris.ix_no_data = ix_no_meris_coverage;
if exist('tca_support.mat','file')==2
    save('tca_support.mat','-append','meris')
else
    save('tca_support.mat','meris')        
end


