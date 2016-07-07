function []=sounding_powerlaw_sens_display()
% function to display the powerlaw sensitivity results.
% results are save in the figures folder.
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
% Bekaert David - university of Leeds
%
% modifications:
% 08/2014   DB  Update the parms_aps file with the values loaded from the
%               earlier processed file
% 03/2016   DB  Remove axis limits and set to default, fix figure size, fix
%               for delay reference height difference 

colorbarflag = 0;
fig_format = 'png';     % 'eps' or 'png'
xlims= [ ];               % empty [] or [start_date end_date], yyyymmdd format 
curdir = pwd;



fontsize=24;
% getting the data from the parms_aps file
sounding_dir = getparm_aps('sounding_dir');
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

% saving of the figures
save_str = [sounding_dir filesep 'Powerlaw' filesep 'figures' filesep];
if exist(save_str,'dir')~=7
    mkdir(save_str)
end

if strcmp(sounding_ifg_dates,'y')
    save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
    n_months =1;
    save_str = [save_str 'SAR_dates'];
    
    
    stamps_processed = getparm_aps('stamps_processed');
    % getting the SAR data information to recompute the InSAR power-law coefficeints
    if strcmp(stamps_processed,'y')   
       fprintf('Stamps processed structure \n')
       ll_matfile = getparm_aps('ll_matfile');
       ps = load(ll_matfile);
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

        keyboard
        fprintf('Check the InSAR dates, this has not been tested \n')
    end
else
    % putting the variables in the right set-up
    start_year = str2num(sounding_start_date(1:4));
    end_year = str2num(sounding_end_date(1:4));
    start_str = sounding_start_date(5:6);
    end_str = sounding_end_date(5:6);

    % the file name to be loaded
    save_name = [sounding_dir filesep 'Powerlaw' filesep 'Powerlaw_sensitivity_' num2str(n_months) 'month_' time_stamp_str 'Hr_' num2str(start_year) start_str '_' num2str(end_year) end_str '.mat'  ];
    save_str = [save_str 'period_'];
end



%loading the data
load([save_name])


% updating the parms_aps file
if strcmp(sounding_ifg_dates,'y')
    alpha_SAR = [alpha_vector_fix(ifgs_ix(:,1)) alpha_vector_fix(ifgs_ix(:,2))];
    h0_SAR = [h0_vector_fix(ifgs_ix(:,1)) h0_vector_fix(ifgs_ix(:,2))];


    % computing the mean between two SAR dates 
    h0_InSAR = mean(h0_SAR,2);
    alpha_InSAR = mean(alpha_SAR,2);

    setparm_aps('powerlaw_h0',h0_InSAR');
    setparm_aps('powerlaw_alpha',alpha_InSAR');
else
     % computing the mean between two SAR dates 
    ix = isnan(h0_vector);
    powerlaw_h0 = mean(h0_vector(~ix));
    powerlaw_alpha = mean(alpha_vector(~ix));

    setparm_aps('powerlaw_h0',powerlaw_h0');
    setparm_aps('powerlaw_alpha',powerlaw_alpha');

end


% getting the mean of the month
date_center = (mean([datenum(date_start_vector,'yyyymmdd') datenum(date_end_vector,'yyyymmdd')],2)) ;

% color variation
color_max = n_months*30*size(time_stamp,1);
color_step = 5;

% generating a colorbar based on the number of soundings per estimate
if isempty(color_max) || color_max==0
	color_max = max(n_soundings_vector);
end
if isempty(color_step)
	color_step = 1;
end

color_max = ceil(color_max/color_step)*color_step;
n_colors = color_max/color_step;


colorcode = jet(n_colors);
color_ix = floor(n_soundings_vector/color_max*n_colors);
color_ix(color_ix==0)=1;

if isempty(xlims)~=1
    xlimits = [datenum(num2str(xlims(1)),'yyyymmdd') datenum(num2str(xlims(2)),'yyyymmdd')];
    window_size = [ 64         254        1493         269];%[30 295 155.6*(max(xlimits)-min(xlimits))./365.25 300];
else
    window_size = [ 64         254        1493         269];%[30 295 155.6*(max(date_center)-min(date_center))./365.25 250];
end

%% Plotting power law coefficient estimated from the 0-hc height range
h_fig = figure('name','Powerlaw from the hc height range','position',window_size);
for k=1:length(alpha_vector)
    if ~isnan(alpha_vector(k))
        h = plot(date_center(k),alpha_vector(k),'ko');
        if colorbarflag==1
            set(h,'Markerfacecolor',colorcode(color_ix(k),:)) ;        
        else
            set(h,'Markerfacecolor','k') ;                    
        end
        hold on
    end
end
xlabel('time of year','fontsize',fontsize)
ylabel('\alpha','fontsize',fontsize)
title('Power law from 0 to hc height range','fontsize',fontsize)
datetick('x','mm/yyyy')
if isempty(xlims)~=1
    xlim(xlimits); 
end
ylim([min(alpha_vector)-(max(alpha_vector)-min(alpha_vector))/10 max(alpha_vector)+(max(alpha_vector)-min(alpha_vector))/10]);

if colorbarflag==1
    cc = colorbar;
    caxis([0 color_max])
    title(cc,'n soundings','fontsize',fontsize)
    colormap(jet(n_colors))
end
set(gca,'fontsize',fontsize)
box on
grid off
set(h_fig,'PaperPositionMode','auto')
if strcmp(fig_format,'eps')
    print(h_fig,'-depsc','-r150',[save_str 'alpha_hc.eps'])
elseif strcmp(fig_format,'png')
    print(h_fig,'-dpng','-r150',[save_str 'alpha_hc.png'])
end

%% Plotting power law coefficient estimated from the full height range
h_fig = figure('name','Powerlaw from the full height range','position',window_size);
for k=1:length(alpha_vector_all)
    if ~isnan(alpha_vector_all(k))
        h = plot(date_center(k),alpha_vector_all(k),'ko');
        if colorbarflag==1
            set(h,'Markerfacecolor',colorcode(color_ix(k),:)) ;        
        else
            set(h,'Markerfacecolor','k');                    
        end
        hold on
    end
end
xlabel('time of year','fontsize',fontsize)
ylabel('\alpha','fontsize',fontsize)
title('Power law from full height range','fontsize',fontsize)
datetick('x','mm/yyyy')
if isempty(xlims)~=1
    xlim(xlimits); 
end
ylim([min(alpha_vector_all)-(max(alpha_vector_all)-min(alpha_vector_all))/10 max(alpha_vector_all)+(max(alpha_vector_all)-min(alpha_vector_all))/10]);
if colorbarflag==1
    cc = colorbar;
    caxis([0 color_max])
    title(cc,'n soundings','fontsize',fontsize)
    colormap(jet(n_colors))
end
set(gca,'fontsize',fontsize)
box on
grid off
set(h_fig,'PaperPositionMode','auto')
if strcmp(fig_format,'eps')
    print(h_fig,'-depsc','-r150',[save_str 'alpha_all.eps'])
elseif strcmp(fig_format,'png')
    print(h_fig,'-dpng','-r150',[save_str 'alpha_all.png'])
end

%% Plotting reference height
h_fig = figure('name','Reference height variation of the power law','position',window_size);
for k=1:length(h0_vector)
    if ~isnan(h0_vector(k))
        h = plot(date_center(k),h0_vector(k),'ko');
        if  colorbarflag==1
            set(h,'Markerfacecolor',colorcode(color_ix(k),:)) ;        
        else
            set(h,'Markerfacecolor','k') ;       
        end
        hold on
    end
end
xlabel('time of year','fontsize',fontsize)
ylabel('h_0 [km]','fontsize',fontsize)
title('Reference height of power law','fontsize',fontsize)
datetick('x','mm/yyyy')
if isempty(xlims)~=1
    xlim(xlimits); 
end
ylim([min(h0_vector)-(max(h0_vector)-min(h0_vector))/10 max(h0_vector)+(max(h0_vector)-min(h0_vector))/10]);

if colorbarflag==1
    cc = colorbar;
    caxis([0 color_max])
    title(cc,'n soundings','fontsize',fontsize)
    colormap(jet(n_colors))
end
set(gca,'fontsize',fontsize)
box on
grid off
set(h_fig,'PaperPositionMode','auto')
if strcmp(fig_format,'eps')
    print(h_fig,'-depsc','-r150',[save_str 'href.eps'])
elseif strcmp(fig_format,'png')
    print(h_fig,'-dpng','-r150',[save_str 'href.png'])
end


%% plotting hte mead delay cure and the mean of the delays

str = '';
while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
    str = input(['Do you want to plot mean delay and net delay? [y/n] \n'],'s');
end
if strcmpi(str,'y')
    fprintf(['Your sounding period is from ' getparm_aps('sounding_start_date') ' till ' getparm_aps('sounding_end_date') '\n'])
    start_date_user = [];
    end_date_user = [];
    while isempty(start_date_user)
        start_date_user = input(['start date yyyymmdd (as a number or NaN for all)?  \n'],'s');
        start_date_user = str2num(start_date_user);
    end
    while isempty(end_date_user)
        end_date_user = input(['end date yyyymmdd (as a number or NaN for all)?  \n'],'s');
        end_date_user = str2num(end_date_user);
    end
    if isnan(start_date_user)
        start_date_user = [];
    else
        start_date_user = num2str(start_date_user);
    end
     if isnan(end_date_user)
        end_date_user = [];
     else
         end_date_user = num2str(end_date_user);
     end   
        
    
    time_stamp = getparm_aps('sounding_time_stamp');
    cd(sounding_dir)

    % making a list of all the sounding files
    [dummy dummy2] = system('echo sounding_list > sounding.list');
    clear dummy dummy2
    for k=1:size(time_stamp,1)
        command_str = ['ls [0-9]???????_' time_stamp(k,:) '.mat >> sounding.list']; 
        [dummy dummy2] = system(command_str);
        clear dummy dummy2
    end

    temp = tdfread('sounding.list');
    [dummy dummy2] = system('rm sounding.list');
    clear dummy dummy2
    date_list_temp = temp.sounding_list(:,[1:8]);
    sounding_list = temp.sounding_list;
    % selecting a date range when requested
    clear ix
       
    abord_flag=0;
    count_runs = 0;
    while abord_flag==0
        if isempty(start_date_user)~=1 && isempty(end_date_user)~=1
            ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date_user,'yyyymmdd') & datenum(date_list_temp,'yyyymmdd')<=datenum(end_date_user,'yyyymmdd'));
        elseif isempty(start_date_user)~=1  && isempty(end_date_user)==1
                ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date_user,'yyyymmdd'));
        elseif isempty(start_date_user)==1 && isempty(end_date_user)~=1
                ix = find(datenum(date_list_temp,'yyyymmdd')<=datenum(end_date_user,'yyyymmdd'));
        end


        if  isempty(start_date_user)~=1 || isempty(end_date_user)~=1
            if isempty(ix)~=1
                date_list_temp = date_list_temp(ix,:);
                sounding_list = sounding_list(ix,:);
            else
                fprintf('No sounding has been acquired in this period \n' )
                continue_flag = 0;
                abord_flag = 1;
            end
            clear temp ix
        end

        % getting the maximum height
        hgt_temp =[];
        for counter_date = 1:size(sounding_list,1)
           % loading of the data
           temp = load(sounding_list(counter_date,:));         
           if isfield(temp,'h_delay');
              hgt_temp = [hgt_temp ;nanmax(temp.h_delay)];
           end
           clear temp
        end
        hgt_temp = sort(hgt_temp);
        max_height_95 = hgt_temp(floor(length(hgt_temp)*0.25));
        max_height_all = max(hgt_temp);
        % check how many we looze for max height based on all heights
        if sum(hgt_temp<max_height_all)>sum(hgt_temp<max_height_95)
            max_height = max_height_95;
        else
            max_height=  max_height_all;
        end
        
        if count_runs==0
           % max_height = 15000;
            height_step = 50;
            max_delay = 200;
            delay_step = 0.25;
            max_sat_delay = [];
            max_sat_netdelay=[];
            scaling_resolution = 1;
        else          
           % max_height = 15000;
            height_step = 50.*scaling_resolution;
            max_delay = 200;
            delay_step = 0.25.*scaling_resolution;
            str = '';
            while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
                str = input(['Are you happy with the color scale or do you want to change manually (y for re-run)? [y/n] \n'],'s');
            end
            abord_flags= 0;
            if strcmpi(str,'n')
                abord_flags = abord_flags+1;            
            else
                h_temp1 = figure('name','color counts'); 
                subplot(2,1,1)
                hist(reshape(delay_matrix_count,[],1))
                title('delays')
                subplot(2,1,2)
                hist(reshape(netdelay_matrix_count,[],1))
                title('netdelays')
                % changing the saturation values
                str = '';
                while isnumeric(str)~=1
                    str = input(['change max saturation for the delay? [] for no change or numeric value \n'],'s');
                    str = str2num(str);
                end
                max_sat_delay = str;
                str= '';
                while isnumeric(str)~=1
                    str = input(['change max saturation for the netdelay? [] for no change or numeric value \n'],'s');
                    str = str2num(str);
                end
                max_sat_netdelay = str;
                
            end
            
            while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
                str = input(['Are you happy with the resolution or do you want to change manually (y for re-run)? [y/n] \n'],'s');
            end
            if strcmpi(str,'n')
                abord_flags = abord_flags+1;
            else
                % changing the saturation values
                str = '';
                while isnumeric(str)~=1
                    str = input(['By what scaling do you want to increase (new = factor*original resolution)? [] for no change or numeric value \n'],'s');
                    str = str2num(str);
                end
                scaling_resolution = str;
                height_step = height_step.*scaling_resolution;
                delay_step = delay_step.*scaling_resolution;
                
                
            end
            
            if abord_flags==2
                abord_flag = 1;
                break
            end
            close(h_delay);
            close(h_netdelay);
            close(h_temp1);
            
        end

 
        
        
        max_height=ceil(max_height./height_step)*height_step;
        max_delay=ceil(max_delay./delay_step)*delay_step;



        nx = [max_delay/delay_step+1];
        ny = [ max_height/height_step+1];
        delay_matrix_count = zeros([ny nx]);
        
        % loop over the different dates and load the data 
        % loading the data for plotting of the delay curves           
        delay_keep=[];
        hgt_keep = [1:max_height]';
        for counter_date = 1:size(sounding_list,1)
            
           % loading of the data
           data_temp = load(sounding_list(counter_date,:));
           
           % looping over the delay 
           if isfield(data_temp,'delay')
%                if ~isempty(delay_keep)
%                    if length(data_temp.delay)==size(delay_keep,1)
%                        
%                        % fix for the reference as some stations go higher
%                        % up than others
%                        ix_reference = find(data_temp.h_delay==max_height);
%                        
%                        delay_keep = [delay_keep data_temp.delay-data_temp.delay(ix_reference)];
%                        sounding_list_keep = [ sounding_list_keep ; sounding_list(counter_date,:)];
% 
%                    end
%                else

                if max(data_temp.h_delay)>=max_height

                    % reference height at the top
                    ix_reference = find(data_temp.h_delay==max_height);

                    % patch from 0 to min sounding height with nan
                    delay_temp = NaN([size(hgt_keep)]);
                    delay_temp(min(data_temp.h_delay):size(delay_temp,1),1) = data_temp.delay(1:ix_reference)-data_temp.delay(ix_reference);
                   
                    delay_keep = [delay_keep delay_temp];   
                    sounding_list_keep = [sounding_list(counter_date,:)];   
                    
                    
                    for m=1:length(hgt_keep)
                        if ~isnan(delay_temp(m))
                          iy_pos = round(hgt_keep(m)./height_step)+1;
                          ix_pos = round(delay_temp(m)*100./delay_step)+1;
                          delay_matrix_count(iy_pos,ix_pos)=delay_matrix_count(iy_pos,ix_pos)+1;
                        end
                   end
                end
                             
           end
           fprintf([num2str(counter_date) ' of ' num2str(size(sounding_list,1)) ' done \n'])
        end 
        if isempty(sounding_list_keep)
            sounding_list_keep=NaN;
        end
        
        fprintf(['Total number of soundings (having data) used in the plot : ' num2str(size(sounding_list_keep,1)) '\n'])

        % computation of the mean delay curve
        mean_delay = mean(delay_keep,2);
        
        % getting the bounding box of the data
        nx_max = max(find(sum(delay_matrix_count,1)~=0));
        ny_max = max(find(sum(delay_matrix_count,2)~=0));
        x_bounds = [0 (nx_max-1)*delay_step]./100;
        y_bounds = [0 (ny_max-1)*height_step]./1000;

        % removing empty space and removing zeros and replace by NaN for plotting
        delay_matrix_count = delay_matrix_count(:,1:nx_max);
        delay_matrix_count = delay_matrix_count(1:ny_max,:);
        delay_matrix_count(delay_matrix_count==0)=NaN;
        
        h_delay = figure('name','Delay curve');
        imagesc(x_bounds,y_bounds,delay_matrix_count)
        axis xy
        if ~isempty(max_sat_delay)
            caxis([0 max_sat_delay]);            
        end
        hold on
        plot(mean_delay,hgt_keep./1000,'k--','linewidth',2) 
        scale_axis =[0 nanmax(nanmax(delay_matrix_count))];
        scale_axis= [floor(scale_axis(1)./1)*1  ceil(scale_axis(2)./1)*1];
        n_steps = (scale_axis(2)-scale_axis(1))./1;
        colormap_figure = (gray(ceil(n_steps*1.15)));
        colormap_figure = flipud(colormap_figure(1:n_steps,:));
        colormap_figure = reshape(repmat(colormap_figure,1,10)',3,[])';
        colormap_figure = [1 1 1 ; colormap_figure];
        colormap(colormap_figure);
        ylabel('h [km]','fontsize',fontsize)
        xlabel('d_{tropo} [m]','fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        
        
        
        %% computation of the net delay
        n_netdelays_max = 400;
        n_soundings = size(delay_keep,2);
        if (n_soundings*n_soundings-n_soundings)/2<=n_netdelays_max
            fprintf(['Make all netdelay differences \n']) 
            ix1_random=repmat(1:n_soundings,n_soundings,1);
            ix1_random=single(ix1_random);
            ix2_random=ix1_random';
            ix_temp=logical(tril(ones(n_soundings)))&~eye(n_soundings);
            ix_random=[ix1_random(ix_temp),ix2_random(ix_temp)]; 
            clear ix1_random ix2_random
        else
            fprintf(['Take random set to make netdelay differences \n']) 
            % take 2 random sets of n_netdelays_max
            new_set_search = 1;
            loop_counter = 0;
            while new_set_search
                ix1_random = ceil(rand(n_netdelays_max*5,1)*n_soundings);
                ix2_random = ceil(rand(n_netdelays_max*5,1)*n_soundings);
                ix_random = [ix1_random ix2_random];
                clear ix1_random ix2_random

                % remove those net delays of the same date
                ix_temp = find(ix_random(:,1)-ix_random(:,2)==0);
                ix_random(ix_temp,:)=[];
                clear ix_temp

                % remove repetition of pair combination
                ix_random = unique(ix_random,'rows');

                % check if enough combinations are left otherwize rerun
                if size(ix_random,1)>=n_netdelays_max
                   ix_random(n_netdelays_max+1:end,:)=[];
                   new_set_search=0;

                end
                loop_counter = loop_counter+1;
                if loop_counter == 50
                    fprintf('Having difficulty determining net delay combinations \n')
                    keyboard
                end
            end
        end

        % computation of the net delay
        netdelay_matrix_temp = [];
        netdelay_matrix_temp = delay_keep(:,ix_random(:,1)) - delay_keep(:,ix_random(:,2));

        
        %%
        netdelay_matrix_count=[];
        max_netdelay = 25;
        delay_step= delay_step./2;

        max_height=ceil(max_height./height_step)*height_step;
        max_netdelay=ceil(max_netdelay./delay_step)*delay_step;

        % initialising of the data matrix
        nx = [2*max_netdelay/delay_step+2];
        ny = [ max_height/height_step+1];
        netdelay_matrix_count = zeros([ny nx]);
        
        
        fprintf(['Total number of netdelays computed for the plot: ' num2str(size(netdelay_matrix_temp,2)) '\n'])
        for counter_date = 1:size(netdelay_matrix_temp,2)
            
           % loading of the data
           data_temp = netdelay_matrix_temp(:,counter_date);
           
           % looping over the delay 
           for m=1:length(data_temp)
               if hgt_keep(m)<=max_height
                   if ~isnan(data_temp(m))
                      iy_pos = round(hgt_keep(m)./height_step)+1;
                      ix_pos = round(data_temp(m)*100./delay_step)+1+nx/2;
                      if ix_pos<=size(netdelay_matrix_count,2) &&  ix_pos>0
                            netdelay_matrix_count(iy_pos,ix_pos)=netdelay_matrix_count(iy_pos,ix_pos)+1;
                      end
                   end
               end                  
           end
           
           if floor(counter_date/50)*50==counter_date
                fprintf([num2str(counter_date) ' of ' num2str(size(netdelay_matrix_temp,2)) ' done \n'])
           end
        end 
        
                
        % getting the bounding box of the data
        nx_max = max(find(sum(netdelay_matrix_count,1)~=0));
        nx_min = min(find(sum(netdelay_matrix_count,1)~=0));
        nx_max = nx/2+max([nx_max-nx/2 nx/2-nx_min]);
        nx_min = nx/2-max([nx_max-nx/2 nx/2-nx_min]);       
        ny_max = max(find(sum(netdelay_matrix_count,2)~=0));
        x_bounds = [(nx_min-nx/2+1)*delay_step (nx_max-nx/2-1)*delay_step]; % cm units
        y_bounds = [0 (ny_max-1)*height_step]./1000;

        if nx_min==0
            nx_min=1;
        end
        % removing empty space and removing zeros and replace by NaN for plotting
        netdelay_matrix_count = netdelay_matrix_count(:,nx_min:nx_max);
        netdelay_matrix_count = netdelay_matrix_count(1:ny_max,:);
        netdelay_matrix_count(netdelay_matrix_count==0)=NaN;
        
        h_netdelay = figure('name','Net-delay curve');
        imagesc(x_bounds,y_bounds,netdelay_matrix_count)
        axis xy
        if ~isempty(max_sat_netdelay)
            caxis([0 max_sat_netdelay]);            
        end
        scale_axis =[0 nanmax(nanmax(netdelay_matrix_count))];
        scale_axis= [floor(scale_axis(1)./1)*1  ceil(scale_axis(2)./1)*1];
        n_steps = (scale_axis(2)-scale_axis(1))./1;
        colormap_figure = (gray(ceil(n_steps*1.2)));
        colormap_figure = flipud(colormap_figure(1:n_steps,:));
        colormap_figure = reshape(repmat(colormap_figure,1,10)',3,[])';
        colormap_figure = [1 1 1 ; colormap_figure];
        colormap(colormap_figure);
        ylabel('h [km]','fontsize',fontsize)
        xlabel('\Delta d_{tropo} [cm]','fontsize',fontsize)
        set(gca,'fontsize',fontsize)

        count_runs = count_runs+1;

    end
    
    cd(curdir)

    
    % plotting the mean delay curve
    str = '';
    while strcmpi(str,'y')~=1 && strcmpi(str,'n')~=1
        str = input(['Do you want to plot the loglog mean delay? [y/n] \n'],'s');
    end
    if strcmpi(str,'y')            
        h0 = getparm_aps('powerlaw_h0');
        % limit to [0 4] km height range as thats topogrpahy
        h_lims = [0 4];
        h_lims_scalled = sort((h0-h_lims));
        
        hgt_keep_new = hgt_keep;
        mean_delay_new = mean_delay;
        mean_delay_new(hgt_keep_new./1000<min(h_lims) |  hgt_keep_new./1000>=max(h_lims))=[];
        hgt_keep_new(hgt_keep_new./1000<min(h_lims) | hgt_keep_new./1000>=max(h_lims))=[];
        h_scalled = (h0-hgt_keep_new./1000);  % scalled heights in km
     
        
        % line through data
        A = [log10(h_scalled) ones(size(h_scalled))];
        coeff = lscov(A,log10(mean_delay_new));
        y_values = [(sort(log10(h_lims_scalled)))' [1;1]]*coeff;
        
        % figure just to get the bounds and the axes for the real figure
        % latter on
        h_temp = figure('position',[1000         392         560         530]);
        plot(log10(h_scalled),log10(mean_delay_new))
        hold on
        plot((sort(log10(h_lims_scalled))),y_values,'r--','linewidth',2)
        
        xlim(sort(log10(h_lims_scalled)))
        xlabel('log(h_0-h [km])','fontsize',fontsize)
        ylabel('log(d_{tropo} [m])','fontsize',fontsize)
        title({'  ','   '},'fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        grid on

        % axis information
        h_ticks = [0.001 1000 2000 3000 4000];
        h_scalled_ticks = log10(h0-h_ticks./1000);  % scalled heights in km
        h_scalled_ticks_label = num2str(round(h_ticks'./1000.*10)./10);
        [h_scalled_ticks,b] = sort(h_scalled_ticks);
        h_scalled_ticks_label = h_scalled_ticks_label(b,:);
        Ax1 = gca;
                

        % start a new figure whit real heights below it and the log height on the top
        figure('position',[1000         392         560         530],'name','Power-law');
        Ax2 = axes('Position',get(Ax1,'Position'));
        plot(log10(h_scalled),log10(mean_delay_new))
        xlim(sort(log10(h_lims_scalled)))
        xlabel('h [km]','fontsize',fontsize)
        set(Ax2,'XAxisLocation','top')
        set(gca,'fontsize',fontsize)
        set(Ax2,'xtick',h_scalled_ticks,'xticklabel',h_scalled_ticks_label)
        
        % put the actual plot on top
        Ax3 = axes('Position',get(Ax1,'Position'));
        hold on
        plot(log10(h_scalled),log10(mean_delay_new),'k-','linewidth',2)
        hold on
        plot((sort(log10(h_lims_scalled))),y_values,'r--','linewidth',2)
        xlim(sort(log10(h_lims_scalled)))
        xlabel('log(h_0-h [km])','fontsize',fontsize)
        ylabel('log(d_{tropo} [m])','fontsize',fontsize)
        title({'  ','   '},'fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        grid on
        box on
        
        % close the dummy figure
        close(h_temp)

    

    
        
    end
    
end


