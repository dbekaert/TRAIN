function []= sounding_powerlaw_sens(hydro,wet)
% function that computes the power law coefficents and performs a
% sensitivity analysis.
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
% By David Bekaert -- University of Leeds 2013
% modifications
% 05/2013   DB:     Compute the power law coefficients for dedicated days.
% 08/2013   DB:     Fix the case when a netdelay cannot be computed.
% 02/2014   DB:     Add non-stamps file support of ifgs-based estimation.
% 03/2014   DB:     Fix the computation of the sensitivity analysis in
%                   absence of soundign for specific SAR dates. 
% 08/2014   DB:     Expand for PS network and SB network for SAR date estimation
% 03/2016   DB:     remove warning on ifgday for non-stamps processing

if nargin<1
    hydro=1;
    wet=1;
end

% getting the data from the parms_aps file
n_months = getparm_aps('sounding_months');
sounding_start_date = getparm_aps('sounding_start_date');
sounding_end_date = getparm_aps('sounding_end_date');
sounding_dir=getparm_aps('sounding_dir');
time_stamp = getparm_aps('sounding_time_stamp');
sounding_ifg_dates = getparm_aps('sounding_ifg_dates');
time_stamp_str = [];
for k=1:size(time_stamp,1)
    if k>1
        time_stamp_str = [time_stamp_str '_' time_stamp(k,:)];
    else
        time_stamp_str = [time_stamp(k,:)];
    end
end
    
% checking if the directory exists
if exist([sounding_dir filesep],'dir')~=7
        error('myApp:argChk', ['The specified filepath of the sounding data does not exist,...  \nAbort,... \n'])
end     
% current directry
curdir = pwd;

% getting the interferogram data information
% loading the data
stamps_processed = getparm_aps('stamps_processed');
if strcmp(stamps_processed,'y')   
   fprintf('Stamps processed structure \n')
   ll_matfile = getparm_aps('ll_matfile');
   ps = load(ll_matfile);
   dates = ps.day;
   lonlat = ps.lonlat;
   
    % getting the parms file list from stamps to see the final ifg list   
    
    % constructing the matrix with master and slave dates
    if strcmp(getparm('small_baseline_flag'),'y')
        % for SB
        ifgs_ix = ps.ifgday_ix;
    
    else
        n_ifg = ps.n_ifg;
        % slightly different for PS.
        date_slave_ix = [1:n_ifg]';
        % the master dates
        date_master_ix = repmat(ps.master_ix,size(date_slave_ix,1),1);
        % ix interferograms
        ifgs_ix = [date_master_ix date_slave_ix];
    end
else
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
    
    for k=1:size(date_master,1)
        [date_master_ix(k,1)] = find(date_master(k,1)==dates);
        [date_slave_ix(k,1)] = find(date_slave(k,1)==dates);
    end
    
    % ix interferograms
    ifgs_ix = [date_master_ix date_slave_ix];
    clear date_master_ix date_slave_ix
end




if strcmp(sounding_ifg_dates,'y')
    % estimate for interferogram dates
    stamps_processed = getparm_aps('stamps_processed');
    if strcmp(stamps_processed,'y')
        ps = load(ll_matfile);
        date_start_vector = datestr(ps.day-15,'yyyymmdd');
        date_end_vector = datestr(ps.day+15,'yyyymmdd');
        
    else

        ifgday_matfile = getparm_aps('ifgday_matfile');
        ifgs_dates = load(ifgday_matfile);
        ifgs_dates = ifgs_dates.ifgday;
        dates = reshape(ifgs_dates,[],1);
        dates = unique(dates);
        dates = datenum(num2str(dates),'yyyymmdd');
        dates = sort(dates);        % dates increasing with time
    
          
        date_start_vector =  datestr(dates-15,'yyyymmdd');
        date_end_vector =  datestr(dates+15,'yyyymmdd');
    end
    
    % save name of the output data
    if hydro==1 && wet==0
        save_name = ['Powerlaw_sensitivity_hydro_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
    elseif hydro==0 && wet==1
        save_name = ['Powerlaw_sensitivity_wet_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
    else
        save_name = ['Powerlaw_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
    end
    if exist([sounding_dir filesep 'Powerlaw'],'dir')~=7
        mkdir([sounding_dir filesep 'Powerlaw']);
    end 
    
    
else
    % estimate is on fixed intervals
    % in case no start or end time is given, get it from the data files
    if isempty(sounding_start_date) || isempty(sounding_end_date)
        
        if exist('sounding.list','file')~=2
            % making a list of all the sounding files
            [dummy dummy2] = system('echo sounding_list > sounding.list');
            clear dummy dummy2
            for k=1:size(time_stamp,1)
                command_str = ['ls [0-9]???????_' time_stamp(k,:) '.mat >> sounding.list']; 
                [dummy dummy2] = system(command_str);
                clear dummy dummy2
            end
        end
        temp = tdfread('sounding.list');
        [dummy dummy2] = system('rm sounding.list');
        clear dummy dummy2
        date_list_temp = temp.sounding_list(:,[1:8]);

        % selecting a date range when requested
        clear ix

        if isempty(start_date)
            sounding_start_date = date_list_temp(1,:);
        end
        if isempty(end_date)
            sounding_end_date = date_list_temp(end,:);    
        end
        clear date_list_temp
    end

    % putting the variables in the right set-up
    start_year = str2num(sounding_start_date(1:4));
    end_year = str2num(sounding_end_date(1:4));
    start_str = sounding_start_date(5:6);
    end_str = sounding_end_date(5:6);
    start_month = str2num(start_str);
    end_month = str2num(end_str);

    % save name of the output data
   if hydro==1 && wet==0
        save_name = ['Powerlaw_sensitivity_hydro_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
    elseif hydro==0 && wet==1
        save_name = ['Powerlaw_sensitivity_wet_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
   else
        save_name = ['Powerlaw_sensitivity_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
    end
    if exist([sounding_dir filesep 'Powerlaw'],'dir')~=7
        mkdir([sounding_dir filesep 'Powerlaw']);
    end 
    
    
    
    if exist([sounding_dir filesep 'Powerlaw'],'dir')~=7
        mkdir([sounding_dir filesep 'Powerlaw']);
    end
    
    % generating the periods
    month_str = ['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'];
    counter = 1;
    % runnign the computation of the powerlaw in batches
    for k=1:end_year-start_year+1
        if k==end_year-start_year+1 & k>1
            for l=1:end_month
                date_start_temp = datenum([num2str(start_year+k-1)  month_str(l,:) '01'],'yyyymmdd');
                date_start_month_vector(counter,:) = datestr(date_start_temp,'yyyymmdd');
                counter = counter+1;
            end
        elseif k==end_year-start_year+1 & k==1
            for l=start_month:end_month
                date_start_temp = datenum([num2str(start_year+k-1)  month_str(l,:) '01'],'yyyymmdd');
                date_start_month_vector(counter,:) = datestr(date_start_temp,'yyyymmdd');
                counter = counter+1;
            end
        else
            for l=1:12
                date_start_temp = datenum([num2str(start_year+k-1)  month_str(l,:) '01'],'yyyymmdd');
                date_start_month_vector(counter,:) = datestr(date_start_temp,'yyyymmdd');
                counter = counter+1;
            end
        end
    end

    n_months_total = size(date_start_month_vector,1);
    for k=1:ceil(n_months_total/n_months)
        k_lower = n_months*(k-1)+1;
        k_upper = n_months*(k)+1;

        if k_upper>n_months_total
            k_upper = n_months_total;
        end
        date_start_vector(k,:)=date_start_month_vector(k_lower,:);
        date_end_vector(k,:) = datestr(datenum(date_start_month_vector(k_upper,:),'yyyymmdd')-1,'yyyymmdd');
    end

end

% Remove those months outside the users request
if ~isempty(sounding_start_date)
    ix_drop = find(datenum(date_start_vector,'yyyymmdd')<datenum(sounding_start_date,'yyyymmdd'));
else
    ix_drop = [];
end
if ~isempty(sounding_end_date)
    ix_drop = [ix_drop ; find(datenum(date_end_vector,'yyyymmdd')>datenum(sounding_end_date,'yyyymmdd'))];
end
ix_drop = unique(ix_drop);
date_end_vector(ix_drop,:)=[];
date_start_vector(ix_drop,:)=[];
    
% Computing the power law coefficients
for k=1:size(date_start_vector,1)
    fprintf(['\n' num2str(k) '/' num2str(size(date_start_vector,1)) ' completed \n']);
	% When having sounding data estimate the powerlaw and b coefficients
    
    [alpha_all,alpha_hc,h_0_threshold,n_soundings] = sounding(date_start_vector(k,:),date_end_vector(k,:),[],hydro,wet);	
	alpha_vector_all(k,1) = alpha_all;               
    alpha_vector(k,1) = alpha_hc;
	n_soundings_vector(k,1)=n_soundings;
	h0_vector(k,1)=h_0_threshold;

end
save([sounding_dir filesep 'Powerlaw' filesep save_name],'alpha_vector_all','alpha_vector','n_soundings_vector','h0_vector','date_start_vector','date_end_vector')

cd(curdir)
if strcmp(sounding_ifg_dates,'y')
    % estimate for interferogram dates
    stamps_processed = getparm_aps('stamps_processed');
    if strcmp(stamps_processed,'y')      
        fprintf('Updating powerlaw parameters with new values. \n')

        % removing NaN by replacing them with the other SAR date estimates
        ix_alpha_fix = find(isnan(alpha_vector(:,1)));
        if ~isempty(ix_alpha_fix)
           if length(ix_alpha_fix)~=length(alpha_vector)
                fprintf([num2str(length(ix_alpha_fix)) ' out of ' num2str(length(alpha_vector)) ' SAR dates where set to the mean alpha due to lack of sounding data. \n']) 
                alpha_vector(isnan(alpha_vector))=nanmean(alpha_vector);
                alpha_vector_fix = alpha_vector;
                save([sounding_dir filesep 'Powerlaw' filesep save_name],'-append','alpha_vector_fix','ix_alpha_fix')

           else
              fprintf('None of the SAR dates where estimate using sounding data \n') 
           end

        end
        
        ix_height_fix = find(isnan(h0_vector(:,1)));
        if ~isempty(ix_height_fix)
           if length(ix_height_fix)~=length(h0_vector)
                fprintf([num2str(length(ix_height_fix)) ' out of ' num2str(length(alpha_vector)) ' SAR dates where set to the mean h0 due to lack of sounding data. \n']) 
                h0_vector(isnan(h0_vector))=nanmean(h0_vector);
                h0_vector_fix = h0_vector;
                save([sounding_dir filesep 'Powerlaw' filesep save_name],'-append','h0_vector_fix','ix_height_fix')
           else
              fprintf('None of the SAR dates where estimate using sounding data \n') 
           end
        end      
        
        
        alpha_SAR = [alpha_vector(ifgs_ix(:,1)) alpha_vector(ifgs_ix(:,2))];
        h0_SAR = [h0_vector(ifgs_ix(:,1)) h0_vector(ifgs_ix(:,2))];
               
        
        % computing the mean between two SAR dates 
        h0_InSAR = mean(h0_SAR,2);
        alpha_InSAR = mean(alpha_SAR,2);
        
        setparm_aps('powerlaw_h0',h0_InSAR');
        setparm_aps('powerlaw_alpha',alpha_InSAR');
    end    
else

    % computing the mean between two SAR dates 
    ix = isnan(h0_vector);
    powerlaw_h0 = mean(h0_vector(~ix));
    powerlaw_alpha = mean(alpha_vector(~ix));

    setparm_aps('powerlaw_h0',powerlaw_h0');
    setparm_aps('powerlaw_alpha',powerlaw_alpha');
end

