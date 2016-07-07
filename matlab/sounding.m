function [alpha_log_all,alpha_log_hc,h_0_threshold,n_soundings] = sounding_new(start_date,end_date,plot_flag,hydro,wet)
% This function computes the refractivity and integrated LOS delays.
% In addition it will estimate the power-law decay cofficient and height at which the relative 
% tropopsheric delays are approximately zero. Both coefficients are given based on the full heifght range
% and the lower height range up to hc. Note that alpha is dimensionless and
% h0 is the height in m units. Sounding data needs to be stored in the 
% correct format! Read below:
% 
% ---- Input ----
% All inputs will be retreived from the aps_parms file. Set the inputs by change 
% setparm_aps('Parameter_name',new_value)
%
% Variables include:
% look_angle     Mean look angle of the data in rad. By default this is
%                assumed to be 21 degees.
% lambda         Radar wavelength in m. By default this is assumed to be
%                0.056 cm (C-band).
% h_0_threshold  Consider only soundings that go up till the h0 threshold in km.
%                Being the height at which the differential delay is small.
%                When '0' is specified it is estimated from the data. Note that
% 		 soundings should go up high enough to cover the lower tropopshere (up to 16 km).
%                By default the threshold is set to be 15 km.
% h_thres_hc     Height range [0 h_thres_hc] over which the power-law decay 
%                coefficient is estimated computed. By default this is set to be 4 km.
% start_date	 Start date of the sounding period, by default [], full period
%                is considered. Specify as a string in 'yyyymmdd' format.
% end_date       End date of the sounding period, by default [], full period 
%                is considered. Specify as a string in 'yyyymmdd' format.
% time_stamp	 Include a time stamp to it. This is a column vector with 
%                strings e.g. ['00';'12'] for 00Z and 12Z. Only those files 
%                ending with this are considered in the computation.
% sounding_dir	 Optional argument giving directly the full path to the soundings.
%                The files on this path should be the YYYYMMDD_HH.mat files. 
% error_promp_flag By default ([] or 'y') this is turn on. Can be usefull to turn of
%                when running in a batch mode. Instead NaN values will be outputed.
%
% Output:
% alpha and h0. The computed sounding delays: refractivity, the (mean) LOS delay [m],
% the (mean) LOS phase delay [rad], and corresponding heights are appended to the existing 
% sounding .mat files. 
%
%
% NOTE on data format:
% sounding data needs to be stored in a folder called sounding_data, 
% within your processing directory. Within this folder, each sounding 
% aquisition needs to be stored in as 8 digit .mat files, e.g. YYYYMMDD_HH.mat 
% format, with the  pressure (hPa), temperature (degree), relative humidity (%)
% and heights (m) as matlab variables P, T, RH and h.
%
% OPTIONAL: 
% - visulize intermediate results (refractivity, delay,
%   mean delay and hc relation with H), by putting plot_flag to 1. By
%   default this is not done.
% - Overwrite all the data to recompute the delays by putting recompute flag to 1
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
% September 2010 --- Bekaert David 
% Modifications:
% 01/11/2011:  DB    Modify from initial version
% 01/12/2011:  DB    Making compatible to other datasets
% 01/11/2012:  DB    Made program more efficient
% 05/11/2012:  DB    Fix bug in loading of the files
% 10/02/2013:  DB    Integrate with other part of toolbox
% 18/03/2013:  DB    Include a time stamp to it.
% 18/03/2013:  DB    Include a start and end time.
% 18/03/2013:  DB    Account for NaN values by removing them for all data columns
% 19/03/2013:  DB    Set the min height range of the delay curve by comparing it to 
%                    the median of the height range
% 19/03/2013:  DB    Include a flag to suppress error messages by lack of sounding data and 
%                    directly pass NaN value through it.
% 25/03/2013:  DB    Compute the power law coefficient from loglog relation
% 03/04/2013:  DB    Include h0 estimation option from the net delays
% 24/04/2013:  DB    Incorporate in the bigger aps toolbox. Change towards
%                    loglog powerlaw computation.
% 10/05/2013:  DB    Loading the default parameters from the parms_aps file.
% 02/10/2013:  DB    Change filename to sounding and clean the syntax
% 29/10/2013:  DB    Allow the look angle to be a specified as a file
% 29/11/2013:  DB    Get stamps processing flag before lookangle  loading
% 06/11/2014:  DB    Add hydrostatic and wet delay options

% --------- VARIABELS ---------- %
fontsize = 20;                      % fontsize of the figures
save_fig =0;                        % when 1 save figures
recompute = 0;                      % Use earlier computed data unless it is missing
if nargin<3 || isempty(plot_flag)
    plot_flag =0;                       % Control plots of the refractivity, delay and net delay.
end
if nargin<4
   hydro=1;
end
if nargin<5
    wet=1;
end
netdelay_color = [];
curdir = pwd;

hydro
wet

%% Setting the defaults were needed
look_angle = getparm_aps('look_angle');
lambda = getparm_aps('lambda');
h_0_threshold= getparm_aps('sounding_h0');
h_thres_hc= getparm_aps('sounding_h_alpha_thres');

if nargin<1 || isempty(start_date)
    start_date= getparm_aps('sounding_start_date');
end 
if nargin<2 || isempty(end_date)
    end_date= getparm_aps('sounding_end_date');
end
time_stamp= getparm_aps('sounding_time_stamp');
sounding_dir= getparm_aps('sounding_dir');
error_promp_flag= getparm_aps('sounding_error_promp');

clear h_0_threshold_default h_thres_hc_default look_angle_default lambda_default

stamps_processed = getparm_aps('stamps_processed');

% checking if the look angle is a file or not
if ischar(look_angle)
	% look angle is specified as a file [rad]
	if strcmp(stamps_processed,'y')
           look_angle = load(look_angle);
	       look_angle = look_angle.la;
           look_angle = mean(look_angle);
	else
	       look_angle = load(look_angle);
	       % use the mean of the look angles
	       look_angle = mean(look_angle);
	end
end

% convert units
theta = look_angle;                       % mean angle of incidence [rad]
h_0_threshold = h_0_threshold*1000;       % put threshold to [m]
h_thres_hc = h_thres_hc*1000;             % put threshold to [m]
clear look_angle

% get height information
if plot_flag==1
    hgt_matfile = getparm_aps('hgt_matfile');
    hgt = load(hgt_matfile);
    max_height = max(hgt.hgt);
    min_height = min(hgt.hgt);
    clear hgt hgt_matfile
end

% Include a while loop. Depending on error message flag, it will break out of the loop and pass
% NaN values as output.
continue_flag =1;
abord_flag = 0;		% when 1 abnormal termination get NaN values instead.
while continue_flag

%% Getting the file list of the sounding data
if isempty(sounding_dir)~=1
	if exist([sounding_dir filesep],'dir')~=7
        error('myApp:argChk', ['The specified filepath of the sounding data does not exist,...  \nAbort,... \n'])
    end     
	cd(sounding_dir)
else
	if exist('sounding_data/','dir')~=7
    		error('myApp:argChk', ['There is no sounding_data directory,...  \nAbort,... \n'])       
	end
	cd  sounding_data
end
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
sounding_list = temp.sounding_list;

% selecting a date range when requested
clear ix

if isempty(start_date)~=1 && isempty(end_date)~=1
	ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date,'yyyymmdd') & datenum(date_list_temp,'yyyymmdd')<=datenum(end_date,'yyyymmdd'));
	save_name = ['Delay_' start_date '_' end_date];
elseif isempty(start_date)~=1  && isempty(end_date)==1
        ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date,'yyyymmdd'));
        save_name = ['Delay_' start_date '_end' ];
elseif isempty(start_date)==1 && isempty(end_date)~=1
        ix = find(datenum(date_list_temp,'yyyymmdd')<=datenum(end_date,'yyyymmdd'));
        save_name = ['Delay_start_' end_date];
else
        save_name = ['Delay'];
end
% checking if the save figures folder exist
if save_fig ==1 & exist('figures','dir')~=7
    mkdir('figures') 
end


if  isempty(start_date)~=1 || isempty(end_date)~=1
	if isempty(ix)~=1
		date_list_temp = date_list_temp(ix,:);
		sounding_list = sounding_list(ix,:);
	else
		fprintf('No sounding has been acquired in this period \n' )
		continue_flag = 0;
                abord_flag = 1;
                break
	end
	clear temp ix
end



%% Computing refractivity
fprintf('Computing the refractivity \n')
ix_skip_sounding = [];			% in case soundings are rejected along the way
warning('off','MATLAB:load:variableNotFound')
for i=1:size(sounding_list,1)
    temp = load(sounding_list(i,:),'h_min_souding');                   
    if isfield(temp,'h_min_souding') && recompute==0
        % storing variables for latter processing
        h_range(i,1) = temp.h_min_souding;          % [m]
        
        if  plot_flag==1
            data = load(sounding_list(i,:),'h_max_souding','Ndelay');
            Ndelay = data.Ndelay;
        else
            data = load(sounding_list(i,:),'h_max_souding');
        end
        h_range(i,2) = data.h_max_souding;          % [m]
        clear data
        skip_sounding = 0;

    else
        load(sounding_list(i,:));                   
        % This includes variables: 
        % - P (pressure in [hPa])
        % - T (temperature in [degrees])
        % - RH (relative humidity in [%]) 
        % - h (altitude in [m])
    
    
        % coping with NaN values in the RH data
        ix1 = find(isnan(RH)==1);
        ix2 = find(isnan(h)==1);
        ix = unique([ix1; ix2]);
        P(ix)=[]; 
        RH(ix)=[];    
        T(ix)=[];        
        h(ix)=[];    
        clear ix ix1 ix2

        % coping with errors in the data when a double recording was made
        ix_repeat = find(diff(sort(h))==0);
        h(ix_repeat) = [];
        P(ix_repeat)=[];
        RH(ix_repeat)=[];        
        T(ix_repeat)=[];
        clear ix_repeat


        % Checking if there is still data left
        if isempty(h)==1 && size(sounding_list,1)==1 && strcmp(error_promp_flag,'y') 
             error('myApp:argChk', ['Datafile contains no numeric data. \n'])
        elseif isempty(h)==1 
            fprintf([sounding_list(i,:) ' sounding skipped as it containes no numeric data\n'])
            skip_sounding = 1;
            ix_skip_sounding = [ix_skip_sounding ; i];
        else
            skip_sounding = 0;
        end

        % when possible compute the refractivity
        if skip_sounding==1
            % case of no data
            h_range(i,1) = NaN;          % [m]
            h_range(i,2) = NaN;          % [m]
        else	
            % saterated pressure es
            Rv = 461.524;                   % [J/K]
            T0 = 273.16;                    % [k]   
            L  = 2.5e6;                     % [J/kg]                Latent heat
            e0 = 6.11;                      % [hPa]             
            T  = T+273.15;                  % [degree] to [K]       Temperature
            es = e0.*exp(L./Rv.*(1./T0-1./T));

            % partial pressure of water vapour e
            RH = RH./100;                   % [%] to [-]            Relative humidity
            e  = RH.*es;                    % [hPa]

            % Refractivety - Hydrostatic term
            k1 = 77.6;                      % [K/hPa]
            Ndelay.Nhydro = k1.*P./T;              % [K/hPa*hPa/K] = [-] 

            % Refractivity - Wet term
            k2_n = 23.3;                    % [K/hPa]
            k3   = 3.75e5;                  % [K^2/hPa]
            Ndelay.Nwet = k2_n.*e./T+k3.*e./T.^2;  % [ppm]

            % refractivety 
            Ndelay.N = Ndelay.Nhydro+Ndelay.Nwet; 

            % save the heights of the delay as well
            Ndelay.h = h;


            h_min_souding = min(h);
            h_max_souding = max(h);

            % saving the data
            save(sounding_list(i,:),'-append','Ndelay','h_min_souding','h_max_souding')

            % storing variables for latter processing
            h_range(i,1) = h_min_souding;          % [m]
            h_range(i,2) = h_max_souding;          % [m]
            clear h_max_souding h_min_souding
            
        end
    end

    if floor(i/10)*10 == i
        fprintf([num2str(i) ' refraction dates completed out of ' num2str(size(sounding_list,1)) ' \n']);
    end
    	
    % Test visualizing refractivity with height
    if  plot_flag==1 && skip_sounding~=1
        if i==1
           figure1 = figure('name','Refractivity with height'); 
           ylabel('Height [km]','fontsize',12)
           xlabel('Refractivity [ppm]','fontsize',12)
           line_colors = hsv(size(sounding_list,1));
           box on
           set(gca,'fontsize',12)
        end
        figure(figure1)
        hold on
        if hydro==1 && wet==0           
            if i==1
                title('Refractivity with height','fontsize',12)
            end
            plot(Ndelay.Nhydro,Ndelay.h./1000,'-','color',line_colors(i,:)) 
        elseif hydro==0 && wet==1
            if i==1
                title('Refractivity with height (hydro)','fontsize',12)
            end
            plot(Ndelay.Nwet,Ndelay.h./1000,'-','color',line_colors(i,:)) 
        else
            if i==1
                title('Refractivity with height (wet)','fontsize',12)
            end
            plot(Ndelay.N,Ndelay.h./1000,'-','color',line_colors(i,:)) 
        end
    end
    clear P h RH T Ndelay N
end
clear Rv T0 L e0 es e k1 k2_n k3 i figure1
% remove those sounding acqusitions that have no data coverage
if isempty(ix_skip_sounding)~=1
	sounding_list(ix_skip_sounding,:)=[];
	h_range(ix_skip_sounding,:)=[];	
end



%% Computing the height range
fprintf('Computing delay \n')

% defining the height range (equal for the whole dataset)
h_max = min(h_range(:,2));
h_max_mean = mean(h_range(:,2));

if h_0_threshold<=0
  % estimate h0 from the data
  h0_estimate_flag = 1;
  h_0_threshold = 0;			% set the height threshold low such all the delays can be used to compute a net delay 
else
  h0_estimate_flag = 0; 
end


% check if the delays goes high enough
if h_max<h_0_threshold
    % the set threshold is lower than the lowest value.
    % setting the h_max to h0, but checking the number of soundings that need to be rejected for this.
    ix = find(h_range(:,2)<h_0_threshold);
    n = [1:size(sounding_list,1)]';
    n(ix)=[];
    h_range(ix,:)=[];
    if isempty(ix)==1
        if strcmp(error_promp_flag,'y') 
            error('myApp:argChk', ['None of the soundings reach higher than ',num2str(h_max/1000),' km. Try setting h_0_threshold lower. \n'])       
        else
            fprintf(['None of the soundings reach higher than ',num2str(h_max/1000),' km. Recommend try to set h_0_threshold lower \n'])
            continue_flag = 0;
            abord_flag = 1;		
            break
        end
    end
    clear ix
    if size(n,1)<3 || isempty(n)==1
        if strcmp(error_promp_flag,'y') 
            error('myApp:argChk', 'Put h_0_threshold to a lower value. \n')
        else
            fprintf('Try to put h_0_threshold to a lower value \n')
            continue_flag = 0; 
            abord_flag = 1;
            break
        end
    end
    if size(n,1)<5
        fprintf(['Only ',num2str(size(n,1)),' soundings are used to compute mean delay \n']);
    end
    sounding_list = sounding_list(n,:);
    h_max = h_0_threshold;
    clear n
elseif h_max>h_0_threshold && h0_estimate_flag==0
    % h0 is lower than the h_max. 
    h_max = h_0_threshold;
elseif h_max_mean-h_max >1000 && h0_estimate_flag==1
    % Make h_max larger possible as the h0 value needs to be estimated
    % this will be at the cost of some soundings but they will be included
    % latter on
    h_max = h_max_mean;
    ix = find(h_range(:,2)<h_max);
    h_range(ix,:)=[];
    sounding_list(ix,:) = [];
end

% including a check for the minimum range by comparing it to the median value
h_min_median = median(h_range(:,1));
h_min = max(h_range(:,1));

% removing soundings that do not start low enough
if h_min-h_min_median > 1000
    h_min = h_min_median;

    ix = find(h_range(:,1)>h_min);
    n = [1:size(sounding_list,1)]';
    n(ix)=[];
    h_range(ix,:)=[];
    if isempty(ix)==1
         error('myApp:argChk', ['This should not occur.'])
    end
    clear ix
    if size(n,1)<3 || isempty(n)==1
        error('myApp:argChk', 'Too few soundings are left for mean delay computation. Check for incomplete data files or differences in sounding height ranges \n.')
    end
    if size(n,1)<5
        fprintf(['Only ',num2str(size(n,1)),' soundings are used to compute mean delay \n']);
    end
    sounding_list = sounding_list(n,:);
    clear n
end

fprintf([num2str(size(sounding_list,1)),' soundings are used to compute mean delay \n']);
h_delay = [h_min:1:h_max]';                                              % set to 1 [m] interval
clear h_range h_min h_max
 

%% computing the integrated refractivity
for i=1:size(sounding_list,1)
    compute_flag=0;
    % checking if it has already been computed. If so use this data at once
    % else compute it
    temp = load(sounding_list(i,:),'h_delay');                   
    if isfield(temp,'h_delay') && recompute==0
        h_delay_temp=temp.h_delay;
        if h_delay(1)==h_delay_temp(1) && h_delay(end)==h_delay_temp(end)
            if hydro==1 && wet==0
                data = load(sounding_list(i,:),'delay_hydro','phase_delay_hydro');
                if ~isfield(data,'delay_hydro')
                   compute_flag=1;
                else
                    delay = data.delay_hydro;
                    phase_delay = data.phase_delay_hydro;
                end
            elseif hydro==0 && wet==1                
                data = load(sounding_list(i,:),'delay_wet','phase_delay_wet');
                if ~isfield(data,'delay_wet')
                    compute_flag=1;
                else
                    delay = data.delay_wet;
                    phase_delay = data.phase_delay_wet;
                end
            else
                data = load(sounding_list(i,:),'delay','phase_delay');
                if ~isfield(data,'delay')
                    compute_flag=1;
                else
                    delay = data.delay;
                    phase_delay = data.phase_delay;
                end
            end
            
            if compute_flag==0
                % store all delays for later mean LOS delay computation
                delay_matrix(1:length(delay),i) = delay;
                % do the same for the phase delay
                phase_delay_matrix(1:length(delay),i) = phase_delay; 
                compute_flag=0;
            end
        else
           % delay needs to be computed
           compute_flag=1;
        end
    else
        compute_flag=1;
    end

    
    if compute_flag==1    
        load(sounding_list(i,:),'Ndelay')
        
        % interpolation to 1 m interval 
        if hydro==1 && wet==0
           N_regular = interp1(Ndelay.h,Ndelay.Nhydro,h_delay,'linear');
        elseif hydro==0 && wet==1
           N_regular = interp1(Ndelay.h,Ndelay.Nwet,h_delay,'linear');
        else
           N_regular = interp1(Ndelay.h,Ndelay.N,h_delay,'linear');
        end
            
            

        % integration of refractivity with height (= delay) and projection on the LOS
        delay = zeros([length(N_regular) 1]);                                   % initialize
        delay = flipud(cumsum(flipud(N_regular))./10^6.*cos(theta));            % 10^6 factor is for correction of ppm, theta for LOS projection
        phase_delay = delay*4*pi./lambda;                                   
        clear scale N_regular 

        % store all delays for later mean LOS delay computation
        delay_matrix(1:length(delay),i) = delay;
        % do the same for the phase delay
        phase_delay_matrix(1:length(delay),i) = phase_delay; 

        % saving the data
        if hydro==1 && wet==0
            delay_hydro = delay;
            phase_delay_hydro = phase_delay;
            save(sounding_list(i,:),'-append','delay_hydro','phase_delay_hydro','h_delay')
        elseif hydro==0 && wet==1
            phase_delay_wet=phase_delay;
            delay_wet = delay;
            save(sounding_list(i,:),'-append','delay_wet','phase_delay_wet','h_delay')
        else
            save(sounding_list(i,:),'-append','delay','phase_delay','h_delay')
        end
        
    end

    % Test visualizing delay with height 
    if plot_flag==1
        if i==1
           figure2 = figure('name','Delay with height'); 
           xlabel('LOS phase delay [rad]','fontsize',fontsize)
           ylabel('Height [km]','fontsize',fontsize)
           set(gca,'fontsize',fontsize)
        end
        figure(figure2)
        hold on
        plot(phase_delay,h_delay./1000,'-','color',[0.7 0.7 0.7])
        box on

        if i==size(sounding_list,1) && h0_estimate_flag==1
            Ax1 = gca;
            ylimits = get(Ax1,'ylim');
            set(Ax1,'ylim',[0 ylimits(2)])
            % adding a second axis 
            box on
            set(gca,'fontsize',fontsize)
            if hydro==1 && wet==0
                title({'LOS hydro delay [m]',' '},'fontsize',fontsize)
                fig_save_name = ['figures' filesep 'Delay_hydro_curves_all.eps'];
            elseif hydro==0 && wet==1
                title({'LOS wet delay [m]',' '},'fontsize',fontsize)
                fig_save_name = ['figures' filesep 'Delay_wet_curves_all.eps'];
            else
                title({'LOS delay [m]',' '},'fontsize',fontsize)
                fig_save_name = ['figures' filesep 'Delay_curves_all.eps'];
            end
            Ax2 = axes('Position',get(Ax1,'Position'),'XAxisLocation','top');
            set(Ax2,'ylim',get(Ax1,'ylim'),'ytick',[],'yticklabel','');
            xlims = get(Ax1,'xlim');
            xtick_new = get(Ax1,'xtick')./xlims(2);
            new_x_label = get(Ax1,'xtick').*lambda./(4*pi);
            xtick_labels = round((new_x_label)*100)./100;
            xtick_labels = num2str(xtick_labels');
            set(Ax2,'xtick',xtick_new,'xticklabel',xtick_labels);
            set(Ax2,'color','none');
            set(Ax2,'fontsize',fontsize)
            box on
            % saving of the figure when requested
            if save_fig==1                
                set(figure2,'PaperPositionMode','auto')
                print(figure2,'-depsc','-r150',fig_save_name)
                clear fig_save_name
            end
            clear figure2 
        end
    end
    clear Ndelay.N Ndelay.h delay 
end
clear i theta

n_soundings = size(delay_matrix,2);


if n_soundings <=3
   continue_flag = 0;
   abord_flag = 1;
   break
end


%% estimating h0 from the delay curves, when required.
% this is based on the net delays and used the std 
% as criteria for setting h0
if h0_estimate_flag==1
    n_netdelays_max = 400;
    if (n_soundings*n_soundings-n_soundings)/2<=n_netdelays_max
       	fprintf(['Estimating h0 from all (', num2str((n_soundings*n_soundings-n_soundings)/2),') possible net delay combinations \n'])
    	ix1_random=repmat(1:n_soundings,n_soundings,1);
        ix1_random=single(ix1_random);
    	ix2_random=ix1_random';
    	ix_temp=logical(tril(ones(n_soundings)))&~eye(n_soundings);
    	ix_random=[ix1_random(ix_temp),ix2_random(ix_temp)]; 
        clear ix1_random ix2_random
    else
        fprintf(['Estimating h0 from (', num2str(n_netdelays_max),') net delay combinations \n'])
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
    netdelays = delay_matrix(:,ix_random(:,1)) - delay_matrix(:,ix_random(:,2));
    netphasedelays = phase_delay_matrix(:,ix_random(:,1)) - phase_delay_matrix(:,ix_random(:,2));

    if plot_flag==1
        fig_netdelay = figure('name','Net delays');
        if isempty(netdelay_color)
            plot(netdelays,h_delay/1000)     
        else
            plot(netdelays,h_delay/1000,'color',netdelay_color)     
        end
        max_spacing = max([abs([min(min(netdelays)) max(max(netdelays))])]);
        xlim([-1.1*max_spacing max_spacing*1.1]);
        xlabel('Net LOS delay [m]','fontsize',fontsize)
        ylabel('Height [km]','fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        
        
        
        fig_netphasedelay = figure('name','Net phase delays');
        if isempty(netdelay_color)
            plot(netphasedelays,h_delay/1000)     
        else
            plot(netphasedelays,h_delay/1000,'color',netdelay_color)     
        end
        max_spacing = max([abs([min(min(netphasedelays)) max(max(netphasedelays))])]);
        xlim([-1.1*max_spacing max_spacing*1.1]);
        xlabel('\Delta\phi_{tropo} [rad]','fontsize',fontsize)
        ylabel('Height [km]','fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        
        
        
        
    end
    
    
    
    
    % computation of the netdelay standard deviation
    netdelays_std = std(netdelays,[],2);
    clear netdelays
    
    % find the height at which the lowest standard deviation becomes larger
    % than the 0.5 cm threshold.
    fprintf(['Estimating h0 from the height at which the std varies more than ' num2str(0.05) ' cm \n'])
    ix = (netdelays_std<=0.0005);
    h_0_threshold = min(h_delay(ix));
    clear ix
    ix = find(h_delay>h_0_threshold);
    
    % make some sanity checks
    if h_0_threshold<6000
       fprintf(['h0 is estimated to be ' num2str(h_0_threshold/1000) ' km. This is low abord. \n'])
       if strcmp(error_promp_flag,'y') 
           error('myApp:argChk', 'Put h_0_threshold to a lower value. \n')
       else
           fprintf('Try setting h_0_threshold manually \n')
           continue_flag = 0; 
           abord_flag = 1;           
           break
       end
    else
        % plotting the estimate h0 on the netdelay curve
        if plot_flag==1
            figure(fig_netdelay)
            hold on
            plot(get(gca,'xlim'),[h_0_threshold/1000 h_0_threshold/1000],'k--')
            
            figure(fig_netphasedelay)
            hold on
            plot(get(gca,'xlim'),[h_0_threshold/1000 h_0_threshold/1000],'k--')
            
            % saving of the figure when requested
            if save_fig==1
                fig_save_name = ['figures' filesep 'NetDelay_curves.eps'];
                set(fig_netdelay,'PaperPositionMode','auto')
                print(fig_netdelay,'-depsc','-r150',fig_save_name)
                clear fig_save_name
                
                fig_save_name = ['figures' filesep 'NetPhaseDelay_curves.eps'];
                set(fig_netphasedelay,'PaperPositionMode','auto')
                print(fig_netphasedelay,'-depsc','-r150',fig_save_name)
                clear fig_save_name
            end
        end
        
        fprintf(['h0 is estimated to be ' num2str(h_0_threshold/1000) ' km. \n'])     

        % removing all the other variables
        h_delay(ix:end)=[];
        delay_matrix(ix:end,:)=[];
        phase_delay_matrix(ix:end,:)=[];
             
        % setting the start of the integration to h_0
        delay_matrix = delay_matrix - repmat(delay_matrix(end,:),length(h_delay),1);
        phase_delay_matrix = phase_delay_matrix - repmat(phase_delay_matrix(end,:),length(h_delay),1);
        
        if plot_flag==1
            figure2 = figure('name','Delay with height'); 
            plot(phase_delay_matrix,h_delay./1000,'-','color',[0.7 0.7 0.7] )    
            xlabel('LOS phase delay [rad]','fontsize',fontsize)
            ylabel('Height [km]','fontsize',fontsize)
            set(gca,'fontsize',fontsize)
            % saving of the figure when requested
            if save_fig==1
                fig_save_name = ['figures' filesep 'Delay_curves.eps'];
                set(figure2,'PaperPositionMode','auto')
                print(figure2,'-depsc','-r150',fig_save_name)
                clear fig_save_name
            end

        end
    end
end 

%% Computing mean delay 
fprintf('Computing mean delay \n')
delay_mean = mean(delay_matrix,2);                                        % given in [m]
phase_delay_mean = mean(phase_delay_matrix,2);
if plot_flag==1
    % plot mean delay on top of the existing figure
    figure(figure2)
    hold on 
    plot(phase_delay_mean,h_delay./1000,'k-','linewidth',2)  
    Ax1 = gca;
    ylimits = get(Ax1,'ylim');
    set(Ax1,'ylim',[0 ylimits(2)])
    % adding a second axis 
    box on
    set(gca,'fontsize',fontsize)
    title({'LOS delay [m]',' '},'fontsize',fontsize)
    Ax2 = axes('Position',get(Ax1,'Position'),'XAxisLocation','top');
    set(Ax2,'ylim',get(Ax1,'ylim'),'ytick',[],'yticklabel','');
    xlims = get(Ax1,'xlim');
    xtick_new = get(Ax1,'xtick')./xlims(2);
    new_x_label = get(Ax1,'xtick').*lambda./(4*pi);
    xtick_labels = round((new_x_label)*100)./100;
    xtick_labels = num2str(xtick_labels');
    set(Ax2,'xtick',xtick_new,'xticklabel',xtick_labels);
    set(Ax2,'color','none');
    set(Ax2,'fontsize',fontsize)
    % saving of the figure when requested
    if save_fig==1
        fig_save_name = ['figures' filesep  'Delay_curves.eps'];
        set(figure2,'PaperPositionMode','auto')
        print(figure2,'-depsc','-r150',fig_save_name)
        clear fig_save_name
    end
    clear figure2

    % an individual figure for the mean delay
    figure4 = figure('name','Mean delay with height'); 
    xlabel('LOS phase delay [rad]','fontsize',fontsize)
    ylabel('Height [km]','fontsize',fontsize)
    set(gca,'fontsize',fontsize)  
    hold on 
    plot(phase_delay_mean,h_delay./1000,'k-','linewidth',2)
    set(gca,'fontsize',fontsize)
end
save([save_name '.mat'],'delay_mean','h_delay','n_soundings','delay_matrix')





%% estimating the power law directly from the log log plot
% computing the loglog powerlaw


h_log = log10(h_0_threshold-h_delay);
phi_delay_log = log10(phase_delay_mean);
delay_log = log10(delay_mean);
A = [h_log ones(size(h_log))];

% based on the full height range
ix_all = find(h_log<=0.001);
A_all = A;
A_all(ix_all,:)=[];
delay_log_temp_all = delay_log;
delay_log_temp_all(ix_all,:)=[];
phase_delay_log_temp_all = phi_delay_log;
phase_delay_log_temp_all(ix_all,:)=[];
coeff_all = inv(A_all'*A_all)*A_all'*delay_log_temp_all;
coeff_all_phase = inv(A_all'*A_all)*A_all'*phase_delay_log_temp_all;
alpha_log_all = coeff_all(1);
K_log_all = exp(coeff_all(2));


% based on the height range till hc_threshold
ix_hc = find(h_log<=0.001 |  h_delay>h_thres_hc);
A_hc = A;
A_hc(ix_hc,:)=[];
delay_log_temp_hc = delay_log;
delay_log_temp_hc(ix_hc,:)=[];
phase_delay_log_temp_hc = phi_delay_log;
phase_delay_log_temp_hc(ix_hc,:)=[];
coeff_hc = inv(A_hc'*A_hc)*A_hc'*delay_log_temp_hc;
coeff_hc_phase = inv(A_hc'*A_hc)*A_hc'*phase_delay_log_temp_hc;
alpha_log_hc = coeff_hc(1);
K_log_hc = exp(coeff_hc(2));



if plot_flag==1
    xlimits = [log10(h_0_threshold-max_height) log10(h_0_threshold+1)];
    % an individual figure for the mean delay
    figure4_1 = figure('name','Log-Log plot of the Mean delay with height');
    ylabel('Log(tropopsheric LOS delay [m])','fontsize',fontsize)
    xlabel('Height [km])','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    hold on
    plot(h_log,delay_log,'k-','linewidth',2)    
    hold on
    plot(A_hc(:,1),A_hc*coeff_hc,'r--','linewidth',2)
    set(gca,'fontsize',fontsize)
    % plotting the powerlaw based on the lower height range
    hold on    
    xlim(xlimits)
    ylimits = get(gca,'ylim');
    legend('Mean delay','Power law',2)
    legend boxoff
    set(gca,'xtick',[get(gca,'xtick') log10(h_0_threshold)],'xticklabel',num2str([get(gca,'xtick') log10(h_0_threshold)]'))

    % updating the lables to be heights
    temp_loc = get(gca,'xtick');
    temp = round((h_0_threshold-10.^temp_loc)/1000*100)/100;
    temp_str = num2str(temp');
    set(gca,'xtick',temp_loc,'xticklabel',temp_str)
    box on
    
    
    if save_fig==1
        fig_save_name = ['figures' filesep 'loglog_mean_delay.eps'];
        set(figure4_1,'PaperPositionMode','auto')
        print(figure4_1,'-depsc','-r150',fig_save_name)
        clear fig_save_name
    end
        
    
    
    
    
    
    
    % an individual figure for the mean delay
    figure4_2 = figure('name','Log-Log plot of the Mean phase delay with height');
    ylabel('Log(\phi_{tropo} [rad])','fontsize',fontsize)
    xlabel('Height [km])','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    hold on
    plot(h_log,phi_delay_log,'k-','linewidth',2)    
    hold on
    plot(A_hc(:,1),A_hc*coeff_hc_phase,'r--','linewidth',2)
    set(gca,'fontsize',fontsize)
    % plotting the powerlaw based on the lower height range
    hold on    
    xlim(xlimits)
    ylimits = get(gca,'ylim');
    legend('Mean delay','Power law',2)
    legend boxoff
    set(gca,'xtick',[get(gca,'xtick') log10(h_0_threshold)],'xticklabel',num2str([get(gca,'xtick') log10(h_0_threshold)]'))

    % updating the lables to be heights
    temp_loc = get(gca,'xtick');
    temp = round((h_0_threshold-10.^temp_loc)/1000*100)/100;
    temp_str = num2str(temp');
    set(gca,'xtick',temp_loc,'xticklabel',temp_str)
    box on
    
    
    if save_fig==1
        fig_save_name = ['figures' filesep 'loglog_mean_phasedelay.eps'];
        set(figure4_2,'PaperPositionMode','auto')
        print(figure4_2,'-depsc','-r150',fig_save_name)
        clear fig_save_name
    end
    
    
    
    
    
    
    
    
    
    
end

save([save_name '.mat'],'-append','alpha_log_all','alpha_log_hc')
clear ix_all A_all delay_log_temp_all coeff_all ix_hc A_hc delay_log_temp_hc coeff_hc delay_log phi_delay_log h_log






%% evaluating the performance of the estimation by plotting the estimated powerlaw wrt to the mean delay from sounding
if plot_flag==1
    % powerlaw estimated from log log plot full height range
    delay_log_all = K_log_all*(h_0_threshold-h_delay).^alpha_log_all;

    % powerlaw delay from log log plot estimated from the 0 till hc height range 
    delay_log_hc = K_log_hc*(h_0_threshold-h_delay).^alpha_log_hc;
    % plot mean delay on top of the existing figure
    hfig_powerlawdelay = figure('name','Fitted power law delay to mean delay curve');
    plot(delay_mean,h_delay./1000,'k-','linewidth',2)
    hold on
    plot(delay_log_hc,h_delay./1000,'r-','linewidth',2)
    hold on
    xlimits = get(gca,'xlim');
    plot([xlimits(1) xlimits(2)],[max_height max_height]./1000,'k--')
    legend('Mean delay','Fitted power law','Max height of the region')
    ylabel('Height [km]','fontsize',fontsize)
    xlabel('LOS delay [m]','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    box on
        
    if save_fig==1
        fig_save_name = ['figures' filesep  'MeanDelay_PowerLawDelay.eps'];
        set(hfig_powerlawdelay,'PaperPositionMode','auto')
        print(hfig_powerlawdelay,'-depsc','-r150',fig_save_name)
        clear fig_save_name
    end
    clear hfig_powerlawdelay
    
end



continue_flag=0;
abord_flag=0;

end

% The while loop was existed because of an error statement.
if abord_flag == 1
	alpha_log_all = NaN;
	alpha_log_hc = NaN;
	n_soundings=NaN;
    h_0_threshold = NaN;
	fprintf('Early termination \n')
else
    h_0_threshold = h_0_threshold/1000;       % put threshold to [km]
	fprintf('DONE \n')
end

cd(curdir)


