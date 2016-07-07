% Script that generates the input files for WPS and WRF model. At the same
% time a file list is generated of the required WRF files. A spin up time
% of 12 Hr is assumed. This can be changed, but its recommended not to
% decrease the spin up time. The GFS data (0.5 deg) is assumed to be on a 6 Hr
% interval. For dates prior to 20061101 GDAS data will be set (1 deg). These are
% based on a 6 Hr interval. 
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
% By Bekaert David  -- University of Leeds
%
% modifications:
% 20/11/2013    DB  Include an input file generation for WPS
% 19/12/2013    DB  Include GDAS files for pre 20061101 dates
% 25/12/2013    DB  Change GDAS to CFSR data, change to surface and pressure data
% 27/12/2013    DB  Include a fix for CFSR as a time stamp earlier is required for WRF
% 06/04/2014    DB  Fixed error in file definition for end of december dates
% 19/09/2014    DB  Add check on dimension for the parent and child grid



%%%--------------------------------%%%
%%%%%% CHANGES REQUIRED BY USER
max_dom = 2;
dx = 25000;
dy = 25000;
parent_grid_ratio =  [1  5];
i_parent_start    =  [1  10];
j_parent_start    =  [1  8];
e_we              = [30  76];
e_sn              = [35  111];
ref_lat           = 18.7;
ref_lon           = -102;
truelat1          = 30.0;
truelat2          = 60.0;
stand_lon         = -102;
map_proj          = 'lambert';
geog_data_path    = '/nfs/a1/software/geog';
geog_data_res     = '3s';


% larger 500 km track:
i_parent_start    =  [1  10];
j_parent_start    =  [1  8];
e_we              = [33  91];
e_sn              = [40  136];

% check for which the WPS software otherwize will crash.
if floor((e_we(2)-1)./parent_grid_ratio(2))~=(e_we(2)-1)./parent_grid_ratio(2) || floor((e_sn(2)-1)./parent_grid_ratio(2))~=(e_sn(2)-1)./parent_grid_ratio(2)
    error(['your grid ratios for child and parent are not good! e_we(child) = parent_ratio * n + 1, where n needs to be integer'])
end

fprintf('Edit aps_wrf_files.m such your domain(s) is(are) set properly \n')
keyboard






%%%%%% NO CHANGES BELOW
%%%--------------------------------%%%


% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed');
ll_matfile = getparm_aps('ll_matfile');
UTC_sat =  getparm_aps('UTC_sat');
t_spin_up = '12:00';                % HH:MM spin up time
data_intervals = 6;                 % data interval in hrs
save_folder = getparm_aps('wrf_datapath');


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
% the weather model spin time:
t_spin_model = t_spin;
t_start_model = t_SAR-t_spin;
% add another time shift to include GFSR data such the full spin up time
% can be completed. Otherzise the first time stamp will not be found.
t_spin = t_spin+data_intervals;
t_start = t_SAR-t_spin;


%% Start and end times
% getting the start date of the WRF model
if t_start_model<0
    n_days_model = floor(t_start_model/24);
    dates_model_wrf = dates-abs(n_days_model);
    time_model_wrf = abs(n_days_model*24)-abs(t_start_model);
else 
    n_days_model = 0;
    dates_model_wrf = dates;
    time_model_wrf = t_start_model;
end
% gettign the start time of the model
time_model_ix_wrf = floor(time_model_wrf/data_intervals);
time_model_wrf = time_model_ix_wrf*data_intervals;
time_model_wrf = num2str(time_model_wrf);
if length(time_model_wrf)==1
   time_model_wrf = ['0' time_model_wrf];
end
% getting the last required time of model data 
time_model_end_wrf = ceil(t_SAR/data_intervals)*data_intervals;
dates_model_end_wrf = dates;
if time_model_end_wrf==24
    time_model_end_wrf=0;
    dates_model_end_wrf = dates_model_end_wrf+1;
end    
time_model_end_wrf = num2str(time_model_end_wrf);
if length(time_model_end_wrf)==1
   time_model_end_wrf = ['0' time_model_end_wrf];
end




% getting the start date for the download of the files.
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


%% Download data for the period:
CFSR_flag =0;
GFS_flag =0;
data_filename_GDAS_prev=[];
if exist([save_folder filesep 'GFS_files.txt'],'file')==2
	delete([save_folder filesep 'GFS_files.txt'])
end
if exist([save_folder filesep 'CFSR_files.txt'],'file')==2
	delete([save_folder filesep 'CFSR_files.txt'])
    delete([save_folder filesep 'CFSR_crop_download.txt'])
end
for k=1:length(dates)
    fprintf(['start_date = ''' datestr(dates_model(k),'yyyy-mm-dd') '_' time_model ':00:00'', end_date = ''' datestr(dates_model_end(k),'yyyy-mm-dd') '_' time_model_end ':00:00'',\n'])

     %% start and end times 
     % for the WRF model
     % getting the model run times in days, hours minutes, seconds
     n_run_days_model_wrf= datenum([datestr(dates_model_end_wrf(k),'yyyymmdd') '_' time_model_end_wrf '0000'],'yyyymmdd_HHMMSS')-datenum([datestr(dates_model_wrf(k),'yyyymmdd') '_' time_model_wrf '0000'],'yyyymmdd_HHMMSS');
     n_run_days_int_model_wrf = floor(n_run_days_model_wrf);
     n_run_hours_model_wrf = (n_run_days_model_wrf-n_run_days_int_model_wrf)*24;
     n_run_hours_int_model_wrf = floor(n_run_hours_model_wrf);
     n_run_minutes_model_wrf = (n_run_hours_model_wrf-n_run_hours_int_model_wrf)*60;
     n_run_minutes_int_model_wrf = floor(n_run_minutes_model_wrf);
     n_run_seconds_model_wrf = round((n_run_minutes_model_wrf-n_run_minutes_int_model_wrf)*60);
     % start time
     start_time_model_wrf = datenum([datestr(dates_model_wrf(k),'yyyymmdd') '_' time_model_wrf '0000'],'yyyymmdd_HHMMSS');
     end_time_model_wrf = datenum([datestr(dates_model_end_wrf(k),'yyyymmdd') '_' time_model_end_wrf '0000'],'yyyymmdd_HHMMSS');
    
    
     % for the downlaod files
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
        
        % GFSR
        GFSR_check = datenum(datestr(data_time,'yyyymmdd'),'yyyymmdd')-datenum('20061101','yyyymmdd'); 
        if GFSR_check<0 && CFSR_flag==0
            fid_files2 = fopen([save_folder filesep 'CFSR_files.txt'],'w');  
            fid_files3 = fopen([save_folder filesep 'CFSR_crop_download.txt'],'w');  

            CFSR_flag=1;
        elseif GFSR_check>=0 && GFS_flag==0
            fid_files = fopen([save_folder filesep 'GFS_files.txt'],'w');
            GFS_flag=1;
        end
        
        % checking if the GFS file exists otherwize use GDAS files
        if GFSR_check<0
            % For GFSR an extra time stamp is needed. 
            % This is not required for GFS data
             skip_flag=0;
             data_types=2;
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

             % zip file contains multiple dates only take those that
             % have not been taken yet.
             if strcmp(data_filename_GDAS,data_filename_GDAS_prev)~=1;
                fprintf(fid_files2,data_filename_GDAS);
                fprintf(fid_files2,data_filename_GDAS2);
                data_filename_GDAS_prev = data_filename_GDAS;
             end
             clear filestr data_filename_GDAS
             
            %% the cropped download script:
            % auto-download cropped for CFSR not full globe download
            startdate_str = [data_time_str(1:4) '-' data_time_str(5:6) '-' data_time_str(7:8) ' ' data_time_str(10:11) ':' data_time_str(12:13)];
            enddate_str = [data_time_str(1:4) '-' data_time_str(5:6) '-' data_time_str(7:8) ' ' data_time_str(10:11) ':' data_time_str(12:13)];

            % 3!7-0.2-1:0.1.1, (relative humidity)
            % 3!7-0.2-1:0.0.0,3!7-4.2-1:0.0.0, (temp)
            %3!7-0.2-1:0.2.3 (v-component of wind)
            %3!7-0.2-1:0.2.2 (u-component of wind)
            %3!7-0.2-1:0.3.1 (pressure reduced to msl)
            %3!7-0.2-1:0.3.5 (geopotential height)

            % defining the crop of the pressure level data
            lon_range = [floor(ref_lon - (dx/1000)*e_we(1)/110)  ceil(ref_lon+(dx/1000)*e_we(1)/110)];
            lon_range(lon_range<-180)=lon_range(lon_range<-180)+360;
            lon_range(lon_range>180)=lon_range(lon_range>180)-360;
            lon_range = sort(lon_range);
            lat_range = [floor(ref_lat - (dy/1000)*e_sn(1)/110)  ceil(ref_lat+(dy/1000)*e_sn(1)/110)];
            lat_range(lat_range<-90)=-90;
            lat_range(lat_range>90)=90;  
            lat_range = sort(lat_range);
            region_string = ['nlat=' num2str(lat_range(2)) ';slat=' num2str(lat_range(1)) ';wlon=' num2str(lon_range(1)) ';elon=' num2str(lon_range(2))  ];   

            % the line used to download the data
            % the actual fetching lines of code
            surface_level_str = ['wget --load-cookies cookiefile --post-data " dsid=ds093.0&rtype=S&rinfo=dsnum=093.0;startdate=' startdate_str ';enddate=' enddate_str ';parameters=3%%217-0.2-1:2.0.192,3%%217-0.2-1:0.3.5,3%%217-0.2-1:0.2.2,3%%217-0.2-1:0.2.3,3%%217-0.2-1:0.1.0,3%%217-0.2-1:0.1.13,3%%217-0.2-1:2.0.0,3%%217-0.2-1:10.2.0,3%%217-0.2-1:0.3.0,3%%217-0.2-1:0.0.0;level=521,522,523,524,107,223,221;product=3;grid_definition=62;ststep=yes;' region_string '" http://rda.ucar.edu/php/dsrqst.php \n'];
            pressure_level_str = ['wget --load-cookies cookiefile --post-data " dsid=ds093.0&rtype=S&rinfo=dsnum=093.0;startdate=' startdate_str ';enddate=' enddate_str ';parameters=3%%217-0.2-1:0.0.0,3%%217-0.2-1:0.1.1,3%%217-0.2-1:0.2.2,3%%217-0.2-1:0.2.3,3%%217-0.2-1:0.3.1,3%%217-0.2-1:0.3.5;level=76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,221,361,362,363,557,562,563,574,577,581,913,914,219;product=3;grid_definition=57;ststep=yes;' region_string '" http://rda.ucar.edu/php/dsrqst.php \n'];

            surface_level_str_test = ['wget -q -O - --post-data " dsid=ds093.0&rtype=S&rinfo=dsnum=093.0;startdate=' startdate_str ';enddate=' enddate_str ';parameters=3%%217-0.2-1:2.0.192,3%%217-0.2-1:0.3.5,3%%217-0.2-1:0.2.2,3%%217-0.2-1:0.2.3,3%%217-0.2-1:0.1.0,3%%217-0.2-1:0.1.13,3%%217-0.2-1:2.0.0,3%%217-0.2-1:10.2.0,3%%217-0.2-1:0.3.0,3%%217-0.2-1:0.0.0;level=521,522,523,524,107,223,221;product=3;grid_definition=62;ststep=yes;' region_string '" http://rda.ucar.edu/php/dsrqst-test.php \n'];
            pressure_level_str_test = ['wget -q -O - --post-data " dsid=ds093.0&rtype=S&rinfo=dsnum=093.0;startdate=' startdate_str ';enddate=' enddate_str ';parameters=3%%217-0.2-1:0.0.0,3%%217-0.2-1:0.1.1,3%%217-0.2-1:0.2.2,3%%217-0.2-1:0.2.3,3%%217-0.2-1:0.3.1,3%%217-0.2-1:0.3.5;level=76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,221,361,362,363,557,562,563,574,577,581,913,914,219;product=3;grid_definition=57;ststep=yes;' region_string '" http://rda.ucar.edu/php/dsrqst-test.php \n'];
            
            if k==1 & kk==1
                user_str = '# Replace <XXX> with the correct information. Do not include < > anymore! \n';
                user_str2 = '# No claims can be made as result of this program -- David Bekaert (University of Leeds) \n';
                coockie_str = 'wget --save-cookies cookiefile --post-data="email=<your_email_address>&passwd=<your_password>&action=login" https://rda.ucar.edu/cgi-bin/login \n';
                fprintf(fid_files3,user_str);
                fprintf(fid_files3,user_str2);
                fprintf(fid_files3,coockie_str);
                clear coockie_str user_str user_str2
            end
            fprintf(fid_files3,surface_level_str);
            fprintf(fid_files3,pressure_level_str);
             
        else
            %% GFS download
            % For GFSR an extra time stamp is needed. 
            % This is not required for GFS data
            if kk>1
                skip_flag = 0;
            else
                skip_flag = 1;
            end
            if skip_flag==0
                data_types=1;
                % GFS file should exist
                data_filename = [datestr(dates(k),'yyyy') filesep 'GFS_Global_0p5deg_' data_time_str '_anl.grib2 \n'];
                fprintf(fid_files,data_filename);
            end
        end
            
    end
    %% closing the writing files
        if k==1:length(dates)
            if skip_flag==0
                if GFS_flag==1
                   fclose(fid_files) 
                end
            end
            if CFSR_flag==1
               fclose(fid_files2) 
               fclose(fid_files3)
            end
        end
    
        
        if max_dom~=0 && skip_flag==0
            SAR_time_min = str2num(UTC_sat(4:5))-str2num(datestr(start_time_model_wrf,'MM'));


            e_we_str = [];
            e_sn_str = [];
            e_vert_str = [];
            dx_vector = dx;
            dy_vector = dy;
            dx_str = [];
            dy_str = [];    
            grid_id_str = [];
            parent_id_str = [];
            i_parent_start_str= [];
            j_parent_start_str= [];
            parent_grid_ratio_str = [];
            input_from_file_str = [];
            frames_per_outfile_str = [];      
            history_interval_str=[]; 
            history_begin_m_str=[];
            nested_str=[];
            specified_str = [];
            start_date_str = [];
            end_date_str = [];
            geog_data_res_str = [];

            % fix the k looping!
            for kk=2:length(parent_grid_ratio)
                dx_vector(kk)= dx_vector(kk-1)/parent_grid_ratio(kk);
                dy_vector(kk)= dy_vector(kk-1)/parent_grid_ratio(kk);
            end
            for kk=1:max_dom
                e_we_str = [e_we_str num2str(e_we(kk)) ', '];
                e_sn_str = [e_sn_str num2str(e_sn(kk)) ', ']; 
                e_vert_str = [e_vert_str num2str(30) ', '];
                dx_str = [dx_str num2str(dx_vector(kk)) ', '];
                dy_str = [dy_str num2str(dy_vector(kk)) ', '];
                grid_id_str = [grid_id_str num2str(kk) ', '];
                parent_id_str = [parent_id_str num2str(kk-1) ', '];
                i_parent_start_str = [i_parent_start_str num2str(i_parent_start(kk)) ', ']; 
                j_parent_start_str = [j_parent_start_str num2str(j_parent_start(kk)) ', ']; 
                parent_grid_ratio_str = [parent_grid_ratio_str num2str(parent_grid_ratio(kk)) ', ']; 
                input_from_file_str = [input_from_file_str '.true., '];
                frames_per_outfile_str = [frames_per_outfile_str num2str(1) ', ']; 
                geog_data_res_str = [geog_data_res_str ' ''' geog_data_res ''', ']; 

                if max_dom==1
                    % when one domain put the output time to correspond to the
                    % SAR time for the parant domain
                    history_begin_m_str = [num2str(SAR_time_min) ', '];
                    history_interval_str = [num2str(60) ', '];
                    nested_str = ['.false.,'];
                    specified_str = ['.true.'];
                    start_date_str = [' ''' datestr(dates_model_wrf(k),'yyyy-mm-dd') '_' time_model_wrf ':00:00'','];
                    end_date_str =  [' ''' datestr(dates_model_end_wrf(k),'yyyy-mm-dd') '_' time_model_end_wrf ':00:00'','];         
                else
                   if kk==1
                       nested_str = ['.false., '];
                       specified_str = ['.true., '];

                   else
                       nested_str = [nested_str '.true., '];
                       specified_str = [specified_str '.false., '];

                   end
                   if kk~=max_dom
                       % when having a nest use only the last time to fix for
                       % the SAR time. Keep parent at the hr timings
                       history_begin_m_str = [history_begin_m_str num2str(0) ', '];
                       % have a 180 min interval for output files generation.
                       history_interval_str = [history_interval_str num2str(180) ', '];
                       % the boundary resets. REset all boudnaries expect for
                       % the smallest domain
                       start_date_str = [start_date_str ' ''' datestr(dates_model_wrf(k),'yyyy-mm-dd') '_' time_model_wrf ':00:00'', '];
                       end_date_str =  [end_date_str ' ''' datestr(dates_model_end_wrf(k),'yyyy-mm-dd') '_' time_model_end_wrf ':00:00'', '];
                   else
                       % the nest with the SAR domain starts output at the SAR
                       % time after the hour on an hourly interval
                       history_begin_m_str = [history_begin_m_str num2str(SAR_time_min) ', '];
                       history_interval_str = [history_interval_str num2str(60) ', '];
                       % the boundary resets. REset all boudnaries expect for
                       % the smallest domain
                       start_date_str = [start_date_str ' ''' datestr(dates_model_wrf(k),'yyyy-mm-dd') '_' time_model_wrf ':00:00'', '];
                       end_date_str =  [end_date_str ' ''' datestr(dates_model_wrf(k),'yyyy-mm-dd') '_' time_model_wrf ':00:00'', '];
                   end

                end
            end

            if exist([save_folder filesep datestr(dates(k),'yyyymmdd')],'dir')~=7
                mkdir([save_folder filesep datestr(dates(k),'yyyymmdd')]);
            end
            fid = fopen([save_folder  filesep datestr(dates(k),'yyyymmdd') filesep 'namelist.input'],'w');


            %%% WRF file        
            fprintf(fid,['\n&time_control\n']) ;
            fprintf(fid,['run_days \t\t = ' num2str(n_run_days_int_model_wrf) ',\n']);
            fprintf(fid,['run_hours \t\t = ' num2str(n_run_hours_int_model_wrf) ',\n']) ;   
            fprintf(fid,['run_minutes \t\t = ' num2str(n_run_minutes_int_model_wrf) ',\n']);
            fprintf(fid,['run_seconds \t\t = ' num2str(n_run_seconds_model_wrf) ',\n']) ;
            fprintf(fid,['start_year \t\t = ' repmat([datestr(start_time_model_wrf,'yyyy') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['start_month \t\t = ' repmat([datestr(start_time_model_wrf,'mm') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['start_day \t\t = ' repmat([datestr(start_time_model_wrf,'dd') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['start_hour \t\t = ' repmat([datestr(start_time_model_wrf,'HH') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['start_minute \t\t = ' repmat([datestr(start_time_model_wrf,'MM') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['start_second \t\t = ' repmat([datestr(start_time_model_wrf,'SS') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_year \t\t = ' repmat([datestr(end_time_model_wrf,'yyyy') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_month \t\t = ' repmat([datestr(end_time_model_wrf,'mm') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_day \t\t = ' repmat([datestr(end_time_model_wrf,'dd') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_hour \t\t = ' repmat([datestr(end_time_model_wrf,'HH') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_minute \t\t = ' repmat([datestr(end_time_model_wrf,'MM') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['end_second \t\t = ' repmat([datestr(end_time_model_wrf,'SS') ', '], 1 , max_dom) '\n']);
            fprintf(fid,['interval_seconds \t = ' num2str(data_intervals*60*60) ',\n']);
            fprintf(fid,['input_from_file \t = ' input_from_file_str '\n']);
            fprintf(fid,['history_interval \t = ' history_interval_str '\n']);
            fprintf(fid,['history_begin_m \t = ' history_begin_m_str '\n']);
            fprintf(fid,['frames_per_outfile \t = ' frames_per_outfile_str '\n']);
            fprintf(fid,['restart \t\t = .false.,\n']);
            fprintf(fid,['restart_interval \t = ' num2str(5000) ',\n']);
            fprintf(fid,['io_form_history \t = ' num2str(2) ',\n']);
            fprintf(fid,['io_form_restart \t = ' num2str(2) ',\n']);
            fprintf(fid,['io_form_input \t\t = ' num2str(2) ',\n']);
            fprintf(fid,['io_form_boundary \t = ' num2str(2) ',\n']);
            fprintf(fid,['debug_level \t\t = ' num2str(0) ',\n']);
            fprintf(fid,['auxhist23_outname \t = "wrfplev_d<domain>_<date>",\n']);
            fprintf(fid,['io_form_auxhist23 \t = ' num2str(2) ',\n']);
            fprintf(fid,['auxhist23_interval \t = ' history_interval_str '\n']);
            fprintf(fid,['frames_per_auxhist23 \t = ' frames_per_outfile_str '\n']);
            fprintf(fid,['auxhist23_begin_m \t = ' history_begin_m_str '\n']);
            fprintf(fid,['/\n\n']);


            fprintf(fid,['\n&domains\n']);
            fprintf(fid,['time_step \t\t = ' num2str(round(mean([dy dx])/1000*6)) ',\n']);
            fprintf(fid,['time_step_fract_num\t = 0,\n']) ;
            fprintf(fid,['time_step_fract_den\t = 1,\n']) ;
            fprintf(fid,['max_dom \t\t = ' num2str(max_dom) '\n']);
            fprintf(fid,['e_we \t\t\t = ' e_we_str '\n']);                          % the grid position for parent and children
            fprintf(fid,['e_sn \t\t\t = ' e_sn_str '\n']);                          % the grid position for parent and children
            fprintf(fid,['e_vert \t\t\t = ' e_vert_str '\n']);
            fprintf(fid,['p_top_requested \t = ' num2str(5000) '\n']);
            if GFSR_check<0
                fprintf(fid,['num_metgrid_levels\t = ' num2str(38) '\n']);         % number of levels for metgrid is different for CFSR
	    else
	    	fprintf(fid,['num_metgrid_levels\t = ' num2str(27) '\n']);         % number of output fileds same as for ERA-I
            end
	    fprintf(fid,['num_metgrid_soil_levels\t = ' num2str(4) '\n']);
            fprintf(fid,['dx \t\t\t = ' dx_str '\n']);
            fprintf(fid,['dy \t\t\t = ' dy_str '\n']);
            fprintf(fid,['grid_id \t\t = ' grid_id_str '\n']);
            fprintf(fid,['parent_id \t\t = ' parent_id_str '\n']);
            fprintf(fid,['i_parent_start \t\t = ' i_parent_start_str '\n']);
            fprintf(fid,['j_parent_start \t\t = ' j_parent_start_str '\n']);
            fprintf(fid,['parent_grid_ratio \t = ' parent_grid_ratio_str '\n']);
            fprintf(fid,['parent_time_step_ratio \t = ' parent_grid_ratio_str '\n']);
            fprintf(fid,['feedback \t\t = ' num2str(1) '\n']);
            fprintf(fid,['smooth_option \t\t = ' num2str(0) '\n']);
            fprintf(fid,['/\n']);


            fprintf(fid,'\n&physics\n');
            fprintf(fid,['mp_physics \t\t = '  repmat('3, ',1,max_dom) '\n']);
            fprintf(fid,['ra_lw_physics \t\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['ra_sw_physics \t\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['radt \t\t\t = ' repmat('30, ',1,max_dom) '\n']);
            fprintf(fid,['sf_sfclay_physics\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['sf_surface_physics\t = '  repmat('2, ',1,max_dom) '\n']);
            fprintf(fid,['bl_pbl_physics\t\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['bldt \t\t\t = '  repmat('0, ',1,max_dom) '\n']);
            fprintf(fid,['cu_physics \t\t = ' repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['cudt \t\t\t = '  repmat('5, ',1,max_dom) '\n']);
            fprintf(fid,['isfflx \t\t\t = 1,\n']);
            fprintf(fid,['ifsnow \t\t\t = 1,\n']);
            fprintf(fid,['icloud \t\t\t = 1,\n']);
            fprintf(fid,['surface_input_source \t = 1,\n']);
            fprintf(fid,['num_soil_layers \t = 4,\n']);
            fprintf(fid,['sf_urban_physics \t = ' repmat('0, ',1,max_dom) '\n']);
            fprintf(fid,['/\n']);

            fprintf(fid,'\n&fdda\n');
            fprintf(fid,['/\n']);

            fprintf(fid,'\n&dynamics\n');
            fprintf(fid,['w_damping  \t\t = 0,\n']);
            fprintf(fid,['diff_opt \t\t = 1,\n']);
            fprintf(fid,['km_opt  \t\t = 4,\n']);
            fprintf(fid,['diff_6th_opt\t\t = '  repmat('0, ',1,max_dom) '\n']);
            fprintf(fid,['diff_6th_factor \t = '  repmat('0.12, ',1,max_dom) '\n']);
            fprintf(fid,['base_temp \t\t = 290,\n']);
            fprintf(fid,['damp_opt \t\t = 0,\n']);
            fprintf(fid,['zdamp\t\t\t = '  repmat('5000., ',1,max_dom) '\n']);
            fprintf(fid,['dampcoef \t\t = '  repmat('0.2, ',1,max_dom) '\n']);
            fprintf(fid,['khdif \t\t\t = '  repmat('0, ',1,max_dom) '\n']);
            fprintf(fid,['kvdif \t\t\t = '  repmat('0, ',1,max_dom) '\n']);
            fprintf(fid,['non_hydrostatic \t = '  repmat('.true., ',1,max_dom) '\n']);
            fprintf(fid,['moist_adv_opt\t\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['scalar_adv_opt\t\t = '  repmat('1, ',1,max_dom) '\n']);
            fprintf(fid,['/\n'])     ;


            fprintf(fid,'\n&bdy_control\n');
            fprintf(fid,['spec_bdy_width  \t = 5,\n']);
            fprintf(fid,['spec_zone  \t\t = 1,\n']);
            fprintf(fid,['relax_zone \t\t = 4,\n']);
            fprintf(fid,['specified \t\t = ' specified_str '\n']);
            fprintf(fid,['nested\t\t\t = ' nested_str '\n']) ;
            fprintf(fid,['/\n']);


            fprintf(fid,'\n&grib2\n');
            fprintf(fid,['/\n']);

            fprintf(fid,'\n&namelist_quilt\n');
            fprintf(fid,['nio_tasks_per_group\t = 0,\n']);
            fprintf(fid,['nio_groups \t\t = 1,\n']);
            fprintf(fid,['/\n']);

            fprintf(fid,'\n&diags\n');
            fprintf(fid,['p_lev_diags \t\t = 1,\n']);
            fprintf(fid,['num_press_levels \t = 37,\n']);
            fprintf(fid,['press_levels \t\t = 100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 75000, 70000, 65000, 60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 22500, 20000, 17500, 15000, 12500, 10000, 7000, 5000, 3000, 2000, 1000, 700, 500, 300, 200, 100,\n']);
            fprintf(fid,['/\n']);

            fclose(fid);






            for ll=1:data_types
                if data_types==2
                   if ll==1
                       fid_wps = fopen([save_folder  filesep datestr(dates(k),'yyyymmdd') filesep 'namelist.wps.pressure'],'w');
                       file_prefix = '''PRESS''';
                       metgrid_prefix = '''SFC'',''PRESS''';

                   elseif ll==2
                       fid_wps = fopen([save_folder  filesep datestr(dates(k),'yyyymmdd') filesep 'namelist.wps.surface'],'w');
                       file_prefix = '''SFC''';
                       metgrid_prefix = '''SFC'',''PRESS''';
                   end
                else
                    fid_wps = fopen([save_folder  filesep datestr(dates(k),'yyyymmdd') filesep 'namelist.wps'],'w');
                    file_prefix = '''FILE''';
                    metgrid_prefix = '''FILE''';
                end
                
                %%% WPS file       
                fprintf(fid_wps,'&share\n');
                fprintf(fid_wps,['wrf_core \t\t\t = ''ARW'',\n']);
                fprintf(fid_wps,['max_dom \t\t\t = ' num2str(max_dom) ',\n']);
                fprintf(fid_wps,['start_date \t\t\t = ' start_date_str '\n']);
                fprintf(fid_wps,['end_date  \t\t\t = ' end_date_str '\n']);
                fprintf(fid_wps,['interval_seconds \t = ' num2str(data_intervals*60*60) ',\n']);
                fprintf(fid_wps,['io_form_geogrid \t = 2,\n']);
                fprintf(fid_wps,['/\n']);

                fprintf(fid_wps,'\n&geogrid\n');
                fprintf(fid_wps,['parent_id \t\t\t = ' parent_id_str '\n']);
                fprintf(fid_wps,['parent_grid_ratio \t = ' parent_grid_ratio_str '\n']);
                fprintf(fid_wps,['i_parent_start \t\t = ' i_parent_start_str '\n']);
                fprintf(fid_wps,['j_parent_start \t\t = ' j_parent_start_str '\n']);
                fprintf(fid_wps,['e_we \t\t\t\t = ' e_we_str '\n']);                          % the grid position for parent and children
                fprintf(fid_wps,['e_sn \t\t\t\t = ' e_sn_str '\n']);                          % the grid position for parent and children
                fprintf(fid_wps,['geog_data_res \t\t = ' geog_data_res_str '\n']);
                fprintf(fid_wps,['dx \t\t\t\t\t = ' num2str(dx) ',\n']);
                fprintf(fid_wps,['dy \t\t\t\t\t = ' num2str(dy) ',\n']);
                fprintf(fid_wps,['map_proj \t\t\t = ''' map_proj ''',\n']);
                fprintf(fid_wps,['ref_lat \t\t\t =  ' num2str(ref_lat) ',\n']);
                fprintf(fid_wps,['ref_lon  \t\t\t =  ' num2str(ref_lon) ',\n']);
                fprintf(fid_wps,['truelat1 \t\t\t = ' num2str(truelat1) ',\n']);
                fprintf(fid_wps,['truelat2 \t\t\t =  ' num2str(truelat2) ',\n']);
                fprintf(fid_wps,['stand_lon \t\t\t =  ' num2str(stand_lon) ',\n']);
                fprintf(fid_wps,['geog_data_path \t\t = ''' geog_data_path  ''',\n']);
                fprintf(fid_wps,['/\n']);

                fprintf(fid_wps,'\n&ungrib\n');
                fprintf(fid_wps,['out_format \t\t\t = ''WPS'',\n']);
                fprintf(fid_wps,['prefix \t\t\t\t = ' file_prefix ',\n']);
                fprintf(fid_wps,['/\n']);

                fprintf(fid_wps,'\n&metgrid\n');
                fprintf(fid_wps,['fg_name \t\t\t = ' metgrid_prefix ',\n']);
                fprintf(fid_wps,['io_form_metgrid \t = 2,\n']);
                fprintf(fid_wps,['/\n']);
                fclose(fid_wps);
            end

        end
end


