function [time_before,time_after, date_before, date_after,f_before,f_after] = aps_weather_model_times(timelist,dates,UTC_sat,model_lag)
% this function gives the time and date for the weather model data before
% and after the SAR acqusition time. It also give the relative weighting to
% the SAR acuisition time of both weather model time-stamps.
% 
% Note in case you have a varying SAR acquisition time, call this script
% iteratively for each UTC time.
% 
% Inputs: 
% timelist is the time at which the weather model data is outputted.
% dates is the julian date of the SAR acquisition(s)
% UTC_sat is the UTC acquisition time of the SAR image.
% 
% KM 02/26/2018  Adding NARR support
if nargin>3
    fprintf('Throwing out dates after model lag time')
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
end
datestructure = 'yyyymmdd';    % assumed date_before(d,:) structure for weather models

% number of SAR acquisitions
n_SAR = length(dates);

time_sampling = 24/(length(timelist)-1);
    
% find two closest times with respect the the 6 hr ERA-I data
time = str2num(UTC_sat(1:2)) + str2num(UTC_sat(4:5))/60;
t_before = floor(time/time_sampling);
t_after = ceil(time/time_sampling);
fprintf(['Satellite pass is ' num2str(time) ' UTC \n'])

% the faction it is closer towards the other date.
f_after = (time - time_sampling*t_before)/(time_sampling*t_after - time_sampling*t_before);
f_after(isnan(f_after))=1;
f_before = 1-f_after;
f_after = repmat(f_after,n_SAR,1);
f_before = repmat(f_before,n_SAR,1);

% the time stamp of the closest two ERA acquisitions
time_before = repmat(num2str(timelist(t_before+1,:)),n_SAR,1);
time_after = repmat(num2str(timelist(t_after+1,:)),n_SAR,1);

% The date for the times after 1800 will change to the next day
clear time

% generating the dates adn times for the before and after SAR acquistion
% weather model files. Accounting for the jump in day.
date_before = datestr(dates,datestructure);
if t_after==length(timelist)-1
    date_after = datestr(dates+1,datestructure);
else
    date_after = datestr(dates,datestructure);
end
