function [] = aps_narr_files(orderflag,model_type)
% Runs and checks which narr data files are needed based on the satellite 
% pass time from the NOAA website and downloads files if 'orderflag' is 1.
% INPUTS: orderflag
% KM February, 2018

% orderflag = 1; % 0: just list data; 1: save data

stdargin = nargin ;
if stdargin<1
    orderflag = 0;
end

model_lag = 8*7;    % days

% Getting the variables from the parms_aps file
workdir = pwd;
stamps_processed = getparm_aps('stamps_processed');
UTC_sat =  getparm_aps('UTC_sat');
narr_datapath = getparm_aps('narr_datapath');
datestructure = 'yyyymmdd';                               % assumed date structure for narr
if isempty(narr_datapath)
    error('please specify narr_datapath')
end

% Loading the data
if strcmp(stamps_processed,'y')
    ll_matfile = getparm_aps('ll_matfile');
    ps = load(ll_matfile);
    dates = ps.day;
    load psver
else
    psver = 2;
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
end

% Make sure all the data exist
end_date = today() - model_lag;
date_fails=dates>end_date;
if sum(date_fails)>0
    fprintf(['The following dates will not have data and are skipped (model lag is ' num2str(model_lag) ')\n'])
    dates_fail = dates(date_fails);
    for k=1:length(dates_fail)
        fprintf([datestr(dates_fail(k),'yyyymmdd') '\n'])
    end
end
% Throw out dates that are later than the lag time of Narr model
dates = dates(~date_fails);

% Get dates
n_dates = length(dates);
% Find two closest times with respect the the 6 hr narr-I data
timelist_NARR = ['0000' ;'0300'; '0600' ; '0900'; '1200' ;'1500'; '1800' ;'2100'; '0000'];
time_sampling = 24/(length(timelist_NARR)-1);
time = str2num(UTC_sat(1:2)) + str2num(UTC_sat(4:5))/60;
t_before = floor(time/time_sampling);
t_after = ceil(time/time_sampling);
fprintf(['Satellite pass is ' num2str(time) ' UTC \n'])
% Ghe fraction it is closer towards the other date.
f_after = (time - time_sampling*t_before)/(time_sampling*t_after - time_sampling*t_before);
f_after(isnan(f_after))=1;
f_before = 1-f_after;
% Ghe time stamp of the closest two narr acquisitions
time1 = num2str(timelist_NARR(t_before+1,:));
time2 = num2str(timelist_NARR(t_after+1,:));
clear time

filelist = [];
clear date date2

% Define two times for each date. If the desired UTC time is after 21:00,
% the second time will come from the subsequent date
for d = 1:n_dates
    date(d,:) = datestr(dates(d),datestructure);
    % Satellite pass is in the evening, next acqusition is next day
    date2(d,:)= date(d,:);
    if t_after==9
        date2(d,:) = datestr(dates(d,:)+1,datestructure);
    end
end

% This is a complete list of all dates and times we need.
date_timelist=([date(:,:), repmat(time1,[n_dates,1]); date2,repmat(time2,[n_dates,1]);]);

% output unique monthly files in case a user needs to find it themselves
datelist=num2str(unique(str2num([date(:,:); date2])));

% This is a list of the unique months we need data for:
monthslist = unique(str2num(datelist(:,1:6))); %This is list of months we need to download

for t = 1:size(monthslist,1)
    file = ['air.' num2str(monthslist(t,:))  '.nc ']; %Format nnn.YYYYMM.nc
    file2 = ['shum.' num2str(monthslist(t,:))  '.nc']; %Format nnnn.YYYYMM.nc
    file3 = ['hgt.' num2str(monthslist(t,:))  '.nc ']; %Format nnn.YYYYMM.nc
    filelist = [filelist ; file;file2;file3];
end

% Saving this information to a file
fid = fopen('NARR_I_files.txt','w');
fid2 = fopen('NARR_I_dates.txt','w');
fprintf(['Required NARR-I files for all interferograms \n'])
for k=1:size(filelist,1)
    fprintf([filelist(k,:) '\n']);
    fprintf(fid,[filelist(k,:) '\n']);
end
fclose(fid);
% Print datelist to console
for k=1:size(datelist,1)
    fprintf(fid2,[datelist(k,:) '\n']);
end
fclose(fid2);

%% Get the NARR data from NOAA
if orderflag==1
    % The URL where the model files exist:
    data_arch = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/';
    % Get the data only if the file doesn't exist 
    lonlat_extracted=0;
    fprintf(['\n Parsing data from ' data_arch ' on NOAA OPeNDAP \n and saving to .mat files in ' narr_datapath ' \n']);

    for ii=1:length(date_timelist) % Loop through each datetime pair
        timestring=strtrim(date_timelist(ii, :));
        month = timestring(1:6);
        day = timestring(1:8);
        dom = str2num(timestring(7:8));
        hhmm= timestring(9:12);

        % Find the index of the desired time within the monthly file 
        timeinteger = find(str2num(timelist_NARR)==str2num(hhmm));
        timeinteger = timeinteger(1); %if it is 00:00, there will be two, so take the first one only.
        nctimeindex = 8*dom -8 + timeinteger; % (8 times/day)

        % Define output names
        outdir=[narr_datapath, filesep, day];
        fout = [outdir, filesep, day, '_', hhmm];
        if ~exist([fout '.mat'],'file')
            if ~exist(outdir,'dir')
                mkdir(outdir);
            end

            % Define file names
            f_air = ['air.', month, '.nc'];
            f_shum = ['shum.', month, '.nc'];
            f_hgt = ['hgt.', month, '.nc'];

            try
                % Parse air temperature file:
                f=[data_arch, f_air];
                if lonlat_extracted==0 %extract metadata only once (doesn't matter which file)
                    lats=ncread(f,'lat'); % grid mapping: lambert_conformal_conic
                    lons=ncread(f,'lon');
                    level=ncread(f,'level');
                    centerlat=ncreadatt(f,'/','centerlat');
                    centerlon=ncreadatt(f,'/','centerlon');
                    latcorners=ncreadatt(f,'/','latcorners');
                    loncorners=ncreadatt(f,'/','loncorners');
                    missing_value=ncreadatt(f, 'air','missing_value');
                    lonlat_extracted=1;
                end
                fprintf(['Saving ' f_air ': ' date_timelist(ii) '\n']);
                air=ncread(f,'air',[1 1 1 nctimeindex],[349 277 29 1]);
                % Format: ncread(filename, 'variable', [startidx startidx startidx startidx], [count count count count]);

                % Parse Specific Humidity file:
                f=[data_arch, f_shum];
                fprintf(['Saving ' f_shum ': ' date_timelist(ii,:) ' \n']);
                shum=ncread(f,'shum',[1 1 1 nctimeindex],[349 277 29 1]);

                % Parse Geopotential Height file
                f=[data_arch, f_hgt];
                fprintf(['Saving ' f_hgt ': ' date_timelist(ii) ' \n']);
                hgt=ncread(f,'hgt',[1 1 1 nctimeindex],[349 277 29 1]);
            catch ME
                fprintf('File not found \n \n')
                fprintf('Matlab failure message: %s\n', ME.message);
            end

            % Save lat, lon, air, shum, and hgt and other metadata to mat file 
            %    as YYYMMDD/YYYMMDD_HHMM.mat        
            save(fout,'air','shum','hgt','lats','lons','level','centerlat','centerlon','latcorners','loncorners','missing_value');
        else
            fprintf([fout ' already exists. Delete before running to get new file. \n'])
        end
    end
else
    fprintf('orderflag is not set to 1. Change if you want to save Narr data.')
end
% Note, to get nc file info: ncdisp(file)