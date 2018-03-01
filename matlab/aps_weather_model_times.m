function [time_before,time_after, date_before, date_after,f_before,f_after,UTC_sat,dates] = aps_weather_model_times(model,ifg_dates,UTC_sat,suppress_output)
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
% modifcations:
% DB    02/2018     supress internal overwrite of timelist


datestructure = 'yyyymmdd';                                 % assumed date_before(d,:) structure for weather models

if nargin<4
    suppress_output='n';
end


% get a list of IFG dates which include the  UTC time for the IFG
[ifg_dates,UTC_sat] = aps_ifg_date_time(ifg_dates,UTC_sat);
data_temp = datestr(ifg_dates(:,1),'yyyymmddHH:MM');
% find the unique set of dates

[data_temp,temp,temp] = unique(data_temp,'rows','sorted');
dates = datenum(data_temp(:,1:8),'yyyymmdd');
UTC_sat = data_temp(:,9:end);
clear temp data_temp ifg_dates 

% getting the dates
n_dates = length(dates);

% retrieve the weather model time spacing
[timelist,n_step] = aps_weather_model_time(model);




%% generate a loop over the number of unique SAR acquistions
time_before=[];
time_after=[];
date_before=[];
date_after=[];
f_before=[];
f_after=[];
for d = 1:n_dates
    time = str2num(UTC_sat(d,1:2)) + str2num(UTC_sat(d,4:5))/60;
    t_before_i = floor(time/n_step);
    t_after_i = ceil(time/n_step);
    if ~strcmpi(suppress_output,'y')
        fprintf(['Satellite pass on ' datestr(dates(d),'yyyymmdd') ' is at ' num2str(time) ' UTC \n'])
    end
    % the faction it is closer towards the other date.
    f_after_i = (time - n_step*t_before_i)/(n_step*t_after_i - n_step*t_before_i);
    f_after_i(isnan(f_after_i))=1;
    f_before_i = 1-f_after_i;

    % the time stamp of the closest two weather model acquisitions
    time_before_i = num2str(timelist(t_before_i+1,:));
    time_after_i = num2str(timelist(t_after_i+1,:));
    % The date for the times after latest model epoch will change to the next day latter on
    clear time

    date_before_i = datestr(dates(d),datestructure);
    
    % Satellite pass is in the evening, next acqusition is next day
    date_after_i= date_before_i;
    if t_after_i==size(timelist,1)-1
        date_after_i = datestr(dates(d)+1,datestructure);
    end
        
 
    % combine the list
    time_before = [time_before;time_before_i];
    time_after = [time_after;time_after_i];
    date_before =[date_before;date_before_i];
    date_after =[date_after;date_after_i];
    f_before =[f_before;f_before_i];
    f_after =[f_after;f_after_i];
end

