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
% By Bekaert David -- University of Leeds
%
% modifications 
% 12/2013 	DB 	Include support for GFSR data
% 12/2013   DB  Add an additional datafile for GFSR such the first time
%               stamp is found. 
% 04/2014   DB  Fixed error in file definition for end of december dates



% Set up the WRF data directory structure
wrf_datapath = [];     % when empty its assumed that data is downloaded to the wrf_path

% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed');
ll_matfile = getparm_aps('ll_matfile');
UTC_sat =  getparm_aps('UTC_sat');
t_spin_up = '12:00';                % HH:MM spin up time
data_intervals = 6;                 % data interval in hrs

if isempty(wrf_datapath)
    wrf_datapath = getparm_aps('wrf_datapath');
end
    
% loading the data
if strcmp(stamps_processed,'y')
   ps = load(ll_matfile);
   dates = ps.day;
   load psver
   fprintf('Stamps processed structure \n')
else
    psver = 2;
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
end
t_SAR = str2num(UTC_sat(1:2))+str2num(UTC_sat(4:5))/60;
t_spin = str2num(t_spin_up(1:2))+str2num(t_spin_up(4:5))/60;
% add another time shift to include GFSR data such the full spin up time
% can be completed. Otherzise the first time stamp will not be found.
t_spin = t_spin+data_intervals;
t_start = t_SAR-t_spin;

% getting the start date of the model
if t_start<0
    n_days = floor(t_start/24);
    dates_model = dates-abs(n_days);
    time_model = abs(n_days*24)-abs(t_start);
else 
    n_days = 0;
    dates_model = dates;
    time_model = t_start;
end
% gettign the start time of the model
time_model_ix = floor(time_model/data_intervals);
time_model = time_model_ix*data_intervals;
time_model = num2str(time_model);
if length(time_model)==1
   time_model = ['0' time_model];
end

% getting the last required time of model data 
time_model_end = ceil(t_SAR/data_intervals)*data_intervals;
dates_model_end = dates;
if time_model_end==24
    time_model_end=0;
    dates_model_end = dates_model_end+1;
end    
time_model_end = num2str(time_model_end);
if length(time_model_end)==1
   time_model_end = ['0' time_model_end];
end


% Download data for the period:
GDAS_flag =0;
GFS_flag =0;
data_filename_GDAS_prev=[];
data_filename_GDAS_prev2=[];
% Download data for the period:
for k=1:length(dates)
    
    % getting the model run times in days, hours minutes, seconds
    n_run_days= datenum([datestr(dates_model_end(k),'yyyymmdd') '_' time_model_end '0000'],'yyyymmdd_HHMMSS')-datenum([datestr(dates_model(k),'yyyymmdd') '_' time_model '0000'],'yyyymmdd_HHMMSS');
    n_run_days_int = floor(n_run_days);
    n_run_hours = (n_run_days-n_run_days_int)*24;
    n_run_hours_int = floor(n_run_hours);
    n_run_minutes = (n_run_hours-n_run_hours_int)*60;
    n_run_minutes_int = floor(n_run_minutes);
    n_run_seconds = round((n_run_minutes-n_run_minutes_int)*60);
     
    % start time
    start_time = datenum([datestr(dates_model(k),'yyyymmdd') '_' time_model '0000'],'yyyymmdd_HHMMSS');
    end_time = datenum([datestr(dates_model_end(k),'yyyymmdd') '_' time_model_end '0000'],'yyyymmdd_HHMMSS');

    n_datafiles = (end_time-start_time)*24/data_intervals+1;
    
    for kk=1:n_datafiles
        data_time = start_time+(kk-1)*data_intervals/24;
        data_time_str = datestr(data_time,'yyyymmdd_HHMM');

        
        if exist([wrf_datapath filesep datestr(dates(k),'yyyymmdd') filesep 'files'],'dir')~=7 
            mkdir([wrf_datapath  filesep datestr(dates(k),'yyyymmdd') filesep 'files']);
        end
  

        % GFSR
        GFSR_check = datenum(datestr(data_time,'yyyymmdd'),'yyyymmdd')-datenum('20061101','yyyymmdd');
        if GFSR_check<0 && GDAS_flag==0
            GDAS_flag=1;
            mkdir([wrf_datapath  filesep 'temp'])
        elseif GFSR_check>=0 && GFS_flag==0
            GFS_flag=1;
        end

        % checking if the GFS file exists otherwize use GDAS files
        if GFSR_check<0
             % For GFSR an extra time stamp is needed. 
             % This is not required for GFS data
             skip_flag=0;
             ix_file_number = ceil(str2num(data_time_str(7:8))/5);
             if ix_file_number==1
                filestr = [data_time_str(1:6) '01-' data_time_str(1:6) '05'];
             elseif ix_file_number==2
                filestr = [data_time_str(1:6) '06-' data_time_str(1:6) '10'];
             elseif ix_file_number==3
                filestr = [data_time_str(1:6) '11-' data_time_str(1:6) '15'];
             elseif ix_file_number==4
                 filestr = [data_time_str(1:6) '16-' data_time_str(1:6) '20'];
             elseif ix_file_number==5
                 filestr = [data_time_str(1:6) '21-' data_time_str(1:6) '25'];
             elseif ix_file_number>=6
                 % get the end day of the month
                 dummy_month = num2str(str2num(data_time_str(5:6))+1);
                 data_time_end_temp = str2num(data_time_str(1:4));
                 
                 if str2num(dummy_month)>12
                     dummy_month = '1';
                     data_time_end_temp = data_time_end_temp+1;
                 end
                 data_time_end_temp = num2str(data_time_end_temp);
                 
                 if length(dummy_month)==1
                     dummy_month = ['0' dummy_month];
                 end
                 end_day_month = datestr(datenum([data_time_end_temp dummy_month '01'],'yyyymmdd')-1,'dd');
                 filestr = [data_time_str(1:6) '26-' data_time_str(1:6) end_day_month];
             end
             data_filename_GDAS = [datestr(dates(k),'yyyy') filesep 'pgbh06.gdas.' filestr '.tar \n'];
             data_filename_GDAS2 = [datestr(dates(k),'yyyy') filesep 'flxf06.gdas.' filestr '.tar \n'];


             if exist([wrf_datapath filesep 'temp'],'file')==2
                 fprintf('There is a file called temp, delete it first otherwize you will loose the downloaded files when performing untar \n')
                 keyboard
             end
                 
             
             % zip file contains multiple dates only take those that
             % have not been taken yet.
             if strcmp(data_filename_GDAS,data_filename_GDAS_prev)~=1;
                data_filename_GDAS_prev = data_filename_GDAS;
                if exist([wrf_datapath  filesep data_filename_GDAS(6:end-3)])
                    command_str = ['tar -xvf ' wrf_datapath  filesep data_filename_GDAS(6:end-2) ' -C ' wrf_datapath  filesep 'temp'];
	                [status,cmdout] = system(command_str);
                    movefile([wrf_datapath  filesep data_filename_GDAS(6:end-3)],[wrf_datapath  filesep 'temp' filesep '.']);
                    command_str2 = ['tar -xvf ' wrf_datapath  filesep data_filename_GDAS2(6:end-2) ' -C ' wrf_datapath  filesep 'temp'];
                    [status,cmdout] = system(command_str2);
                    if exist([wrf_datapath  filesep data_filename_GDAS2(6:end-3)],'file')==2
			     movefile([wrf_datapath  filesep data_filename_GDAS2(6:end-3)],[wrf_datapath  filesep 'temp' filesep '.']);
		    end
                end
            end
       	    if strcmp(data_filename_GDAS2,data_filename_GDAS_prev2)~=1;
                data_filename_GDAS_prev2 = data_filename_GDAS2;
                if exist([wrf_datapath  filesep data_filename_GDAS2(6:end-3)])
                        command_str2 = ['tar -xvf ' wrf_datapath  filesep data_filename_GDAS2(6:end-2) ' -C ' wrf_datapath  filesep 'temp'];      
                       [status,cmdout] = system(command_str2);
			if exist([wrf_datapath  filesep data_filename_GDAS2(6:end-3)],'file')==2
	                        movefile([wrf_datapath  filesep data_filename_GDAS2(6:end-3)],[wrf_datapath  filesep 'temp' filesep '.']);
			end
                end
            end
            data_time_str_new = datestr(datenum(data_time_str(1:end-2),'yyyymmdd_HH'),'yyyymmddHH');
            if exist([wrf_datapath filesep 'temp' filesep 'pgbh06.gdas.' data_time_str_new '.grb2'])
                if exist([wrf_datapath  filesep datestr(dates(k),'yyyymmdd') filesep 'files' filesep 'pressure'],'dir')~=7
                    mkdir([wrf_datapath  filesep datestr(dates(k),'yyyymmdd') filesep 'files' filesep 'pressure']);
                end
		if exist([wrf_datapath filesep 'temp' filesep 'pgbh06.gdas.' data_time_str_new '.grb2'],'file')==2
	                movefile([wrf_datapath filesep 'temp' filesep 'pgbh06.gdas.' data_time_str_new '.grb2'],[wrf_datapath filesep  datestr(dates(k),'yyyymmdd')  filesep 'files' filesep 'pressure' filesep 'pgbh06.gdas.' data_time_str_new '.grb2'])
		end
            end
            if exist([wrf_datapath filesep 'temp' filesep 'flxf06.gdas.' data_time_str_new '.grb2'])                
                if exist([wrf_datapath  filesep datestr(dates(k),'yyyymmdd') filesep 'files' filesep 'surface'],'dir')~=7
                      mkdir([wrf_datapath  filesep datestr(dates(k),'yyyymmdd') filesep 'files' filesep 'surface']);                 
                end
		if exist([wrf_datapath filesep 'temp' filesep 'flxf06.gdas.' data_time_str_new '.grb2'],'file')==2
	            movefile([wrf_datapath filesep 'temp' filesep 'flxf06.gdas.' data_time_str_new '.grb2'],[wrf_datapath filesep  datestr(dates(k),'yyyymmdd')  filesep 'files' filesep 'surface' filesep 'flxf06.gdas.' data_time_str_new '.grb2'])
		end
            end        

    
             clear filestr data_filename_GDAS
        else
            % GFS file should exist 
            % For GFSR an extra time stamp is needed. 
            % This is not required for GFS data
            if kk>1
                skip_flag = 0;
            else
                skip_flag = 1;
            end
            if skip_flag==0
                if exist([wrf_datapath filesep 'GFS_Global_0p5deg_' data_time_str '_anl.grib2'])     
            		movefile([wrf_datapath filesep 'GFS_Global_0p5deg_' data_time_str '_anl.grib2'],[wrf_datapath filesep  datestr(dates(k),'yyyymmdd')  filesep 'files' filesep 'GFS_Global_0p5deg_' data_time_str '_anl.grib2'])
                end	
            end
        end
     end

end
    
