function [time_before,time_after, date_before, date_after,f_before,f_after] = aps_weather_model_times(timelist,dates,UTC_sat)
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


datestructure = 'yyyymmdd';                                 % assumed date_before(d,:) structure for weather models

% number of SAR acquisitions
n_SAR = length(dates);

% find two closest times with respect the the 6 hr ERA-I data
timelist = ['0000' ; '0600' ; '1200' ; '1800' ; '0000'];
time = str2num(UTC_sat(1:2)) + str2num(UTC_sat(4:5))/60;
t_before = floor(time/6);
t_after = ceil(time/6);
fprintf(['Satellite pass is ' num2str(time) ' UTC \n'])

% the faction it is closer towards the other date.
f_after = (time - 6*t_before)/(6*t_after - 6*t_before);
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
if t_after==4
    date_after = datestr(dates+1,datestructure);
else
    date_after = datestr(dates,datestructure);
end
