function [] = sounding_spectrometer_sens(k1,k2,k3,h_max)
% Program which extract the surface temperature from the sounding list and 
% computes the PI conversion factor for spectrometer estimation of the tropopsheric delay.
%
% OPTIONAL INPUTS
% k1            Optional argument to specify k1 of the soundign eq in [K/Pa]
% k2            Optional argument to specify k2 of the soundign eq in [K/Pa]
% k3            Optional argument to specify k3 of the soundign eq in [K^2/Pa]
% h_max         optional argument that gives they height for which below
%               the mean temperature is computed for the atmopshere, used  
%               in the computation of the scale height 
%
% INPUTS loaded from the parms_aps file
% start_date	Start date of the sounding period, by default [], full period
%               is considered.
% end_date      End date of the sounding period, by default [], full period 
%               is considered.
% time_stamp	Include a time stamp to it. This is a column vector with 
%               strings e.g. ['00';'12'] for 00Z and 12Z. Only those files 
%               ending with this are considered in the computation.

% sounding_dir	Optional argument giving directly the full path to the soundings.
%               The files on this path should be the YYYYMMDD_HH.mat files. 
% error_promp_flag By default ([] or 1) this is turn on. Can be usefull to turn of
%               when running in a batch mode. Instead NaN values will be outputed.
%
% Output:
% Creates a .mat files with the PI conversion factor and surface temperature associated
% with different dates. 
%
%
% NOTE on data format:
% sounding data needs to be stored in a folder called sounding_data, 
% within your processing directory. Within this folder, each sounding 
% aquisition needs to be stored in as 8 digit .mat files, e.g. YYYYMMDD_HH.mat 
% format, with the  pressure (hPa), temperature (degree), relative humidity (%)
% and heights (m) as matlab variables P, T, RH and h.
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
% April 2013 --- Bekaert David 
% Modifications:
% 12/04/2013    Bekaert David   Include compuation for PI factor from
%                               Richard Walters shell script.
% 21/10/2013    Bekaert David   Include in APS toolbox, retrieve parms
%                               from parms_aps.
% 19/03/2014    Bekaert David   Suppress command line output, same setup of
%                               variable saving as for power-law
% 14/05/2014	Bekaert David	rename to spectrometer function name

% --------- VARIABELS ---------- %
fontsize = 15;                      % fontsize of the figures
save_fig = 1;

%% Setting the defaults were needed
% setting the constants 
if nargin<1 || isempty(k1)
    k1 = 0.776;                     % [K/Pa]
end 
if nargin<2 || isempty(k2)
    k2 = 0.704;                     % [K/Pa]
end
if nargin<3 || isempty(k3)
    k3   = 3739;                    % [K^2/Pa]
end
if nargin<4 || isempty(h_max)
   h_max = 20000;                   % mean temp of heighs below this maximum 
                                    % height are used to compute the mean
                                    % temperature for the scale height
end

curdir=pwd;

% getting parameters from parms_aps file
sounding_dir = getparm_aps('sounding_dir');
error_promp_flag = getparm_aps('sounding_error_promp');

% getting teh time stamp needed
time_stamp = getparm_aps('sounding_time_stamp');
time_stamp_str = [];
for k=1:size(time_stamp,1)
    if k>1
        time_stamp_str = [time_stamp_str '_' time_stamp(k,:)];
    else
        time_stamp_str = [time_stamp(k,:)];
    end
end

% check if it needs to be computed on SAR dates or not
sounding_ifg_dates = getparm_aps('sounding_ifg_dates');
if strcmp(sounding_ifg_dates,'y')
    sar_date_sounding =1;
    
    % estimate for interferogram dates
    % It used a window arounf the data, but it wil first check if the date
    % itself has data and use that value. Incase there is no data for the
    % SAR date it will use an average 
    
    % loading the data    
    stamps_processed = getparm_aps('stamps_processed');
    if strcmp(stamps_processed,'y')
       ll_matfile = getparm_aps('ll_matfile',1);
       ps = load(ll_matfile);
       ifgs_dates = ps.day;
       fprintf('Stamps processed structure \n')
    else
        ifgday_matfile = getparm_aps('ifgday_matfile',1);
        ifgs_dates = load(ifgday_matfile);
        ifgs_dates = ifgs_dates.ifgday;
        ifgs_dates = reshape(ifgs_dates,[],1);
        ifgs_dates = unique(ifgs_dates);
    end

    date_start_vector =  datestr(ifgs_dates-15,'yyyymmdd');
    date_end_vector =  datestr(ifgs_dates+15,'yyyymmdd');
    
    % save name of the output data
    date_start_temp = getparm_aps('sounding_start_date');
    date_end_temp = getparm_aps('sounding_end_date');
    start_year_str = date_start_temp(1:4);
    end_year_str = date_end_temp(1:4);
    start_str = date_start_temp(5:6);
    end_str = date_end_temp(5:6);
    
    save_name_final = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat']

    
    if exist([sounding_dir filesep 'Spectrometer'],'dir')~=7
        mkdir([sounding_dir filesep 'Spectrometer']);
    end 
    
    
else
    % this is regular sensitivity on a wider range of data
    sar_date_sounding =0;
    date_start_vector = getparm_aps('sounding_start_date');
    date_end_vector = getparm_aps('sounding_end_date');
   
    start_year_str = date_start_vector(1:4);
    end_year_str = date_end_vector(1:4);
    start_str = date_start_vector(5:6);
    end_str = date_end_vector(5:6);
    
    save_name_final = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat'];
end


% Constants
R = 8.3144621;      % [J K-1 mol-1]
Mw = 0.01801528;    % [kg/mol]
Md = 0.02897;       % moleculair mass of air [kg/mol]
pw = 1000;          % density of water [kg/m3]
g = 9.80666;        % gravitational acceleration [m/s^2]

if save_fig ==1 & exist([sounding_dir filesep 'Spectrometer' filesep 'figures'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer' filesep 'figures']) 
end
if save_fig ==1 & exist([sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor']) 
end
if save_fig ==1 && exist([sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'H_scaling_factor'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor']) 
end



% looping over the data. In case SAR sensitivity is needed no average is
% pased:
for date_counter=1:size(date_start_vector,1)

    start_date = date_start_vector(date_counter,:);
    end_date = date_end_vector(date_counter,:);


    start_year_str = start_date(1:4);
    end_year_str = end_date(1:4);
    start_str = start_date(5:6);
    end_str = end_date(5:6);

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
            command_str = ['echo sounding_list > sounding.list']; 
            [a,b] =system(command_str);
            clear command_str
            for k=1:size(time_stamp,1)
                command_str = ['ls [0-9]???????_' time_stamp(k,:) '.mat >> sounding.list']; 
                [a,b] =system(command_str);
            end
        end
        temp = tdfread('sounding.list');
        command_str = ['rm sounding.list']; 
        [a,b] =system(command_str);
        clear command_str
        
        date_list_temp = temp.sounding_list(:,[1:8]);
        sounding_list = temp.sounding_list;

        % selecting a date range when requested
        clear ix
        if isempty(start_date)~=1 && isempty(end_date)~=1
            ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date,'yyyymmdd') & datenum(date_list_temp,'yyyymmdd')<=datenum(end_date,'yyyymmdd'));
        elseif isempty(start_date)~=1  && isempty(end_date)==1
                ix = find(datenum(date_list_temp,'yyyymmdd')>=datenum(start_date,'yyyymmdd'));
        elseif isempty(start_date)==1 && isempty(end_date)~=1
                ix = find(datenum(date_list_temp,'yyyymmdd')<=datenum(end_date,'yyyymmdd'));
        else
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
        fprintf('Loading the data \n')
        ix_skip_sounding = [];			% in case soundings are rejected along the way

        % initialisation of the variables
        hs_vector = NaN([size(sounding_list,1) 1]);
        date_vector_J = NaN([size(sounding_list,1) 1]);
        Ts_vector = NaN([size(sounding_list,1) 1]);
        Tm_vector = NaN([size(sounding_list,1) 1]);
        h_scaling_vector = NaN([size(sounding_list,1) 1]);
        T_mean_vector = NaN([size(sounding_list,1) 1]);
        Pi_factor_vector = NaN([size(sounding_list,1) 1]);

        for i=1:size(sounding_list,1)
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
            if isempty(h)==1 && size(sounding_list,1)==1 && error_promp_flag==1 
                 error('myApp:argChk', ['Datafile contains no numeric data. \n'])
            elseif isempty(h)==1 
                fprintf([sounding_list(i,:) ' sounding skipped as it containes no numeric data\n'])
                skip_sounding = 1;
                ix_skip_sounding = [ix_skip_sounding ; i];
            else
                skip_sounding = 0;
            end

            % when possible compute the PI factor
            if skip_sounding==1
                % case of no data
                hs = NaN;          % surface elevation [m]
                Ts = NaN;          % surface temperature [deg]
                Tm = NaN;
                PI = NaN;

            else	
                % getting the surface elevation and surface temperature
                hs = h(1);          % surface elevation [m]
                Ts = T(1);          % surface temperature [deg]

                % computation of the PI factor
                Tm = (70.2 + 0.72*(Ts + 273.15)) -273.15;
                PI = 1e-6* pw*(R/Mw) * (k3/(Tm+ 273.15) + (k2 -(k1*(Mw/ Md))));
            end
            PI_factor.PI = PI;
            PI_factor.Tm = Tm;
            PI_factor.Ts = Ts;
            PI_factor.hs = hs;

            % computation of the scaling height
            T_mean  = mean(T(h<=h_max));
            H_scaling_factor = R*(T_mean+273.15)/(Md*g);
            H_scaling.H_scaling = H_scaling_factor;
            H_scaling.T_mean = T_mean;


            % saving the data
            save(sounding_list(i,:),'-append','PI_factor','H_scaling')   

            % Generating a vector for the output
            Pi_factor_vector(i)=PI;
            date_vector_J(i)=datenum(date_list_temp(i,:),'yyyymmdd');
            hs_vector(i)=hs;
            Ts_vector(i)=Ts;
            Tm_vector(i)=Tm;
            h_scaling_vector(i) = H_scaling_factor;
            T_mean_vector(i) = T_mean;


            clear P h RH T PI Tm TS hs PI_factor H_scaling_factor T_mean

            if floor(i/10)==i/10
                fprintf([num2str(i) ' out of ' num2str(size(sounding_list,1)) ' done \n'])
            end
        end
        clear Rv T0 L e0 es e i figure1

        % remove those sounding acqusitions that have no data coverage
        if isempty(ix_skip_sounding)~=1
             sounding_list(ix_skip_sounding,:)=[];
        end

        continue_flag=0;
        abord_flag=0;
    end

    % The while loop was existed because of an error statement.
    if abord_flag == 1 && sar_date_sounding ==0
        Pi_factor_vector=NaN;
        date_vector_J=NaN;
        hs_vector=NaN;
        Ts_vector=NaN;
        Tm_vector=NaN;
        h_scaling_vector = NaN;
        T_mean_vector = NaN;
        fprintf('Early termination \n')
    elseif abord_flag == 1 && sar_date_sounding ==1
        spectrometer_scaleheight(date_counter)=NaN;
        spectrometer_PIconversion(date_counter)=NaN;
    else 
        % checking for outliers
        ix_outlier = [];
        % surface elevations. It might occur that the sounding ballon only
        % start acquiring at a latter stage
        ix_nan = isnan(hs_vector);
        hs_median = median(hs_vector(~ix_nan));
        % search for those sounding where the elevation differs by 20 m from
        % the median elevation
        ix_hs = find(abs(hs_vector-hs_median)>=20); 
        fprintf('\n Search for soundings that did not start acquiring at the surface \n')

        if isempty(ix_hs)
            fprintf('All soundings have their surface height within +-20m of the median elevation. \n')
        else
            fprintf([num2str(length(ix_hs)) ' soundings were found to deviate more than +-20m of the median elevation.\n'])
            ix_outlier = [ix_outlier ; ix_hs];
        end
        % Mean temperatures. It might occur that the sounding ballon only
        % start acquiring at a latter stage
        Ts_mean = mean(Ts_vector(~ix_nan));
        % search for those sounding where the temperature differs by 15 deg from
        % the median elevation
        ix_Ts = find(abs(Ts_vector-Ts_mean)>15);      
        if isempty(ix_Ts)
            fprintf('All soundings have their mean temperatue within +-15deg of the mean temperature of all soundings. \n')
        else
            fprintf([num2str(length(ix_Ts)) ' soundings were found to have their mean temperatue to deviate more than +-15deg of the mean temperature of all soundings.\n'])
            ix_outlier = [ix_outlier ; ix_Ts];
        end
        % remove potential duplicates
        ix_outlier=unique(ix_outlier);

        if save_fig==1 && sar_date_sounding ==0
            if isempty(ix_outlier)~=1
                hfig = figure('name','Pi-factor variation with potential outliers');
                plot(date_vector_J,Pi_factor_vector,'k.')
                datetick('x','mmm/yyyy')
                ylabel('PI Factor','fontsize',fontsize)
                title(['PI-factor for surface temperature at ' num2str(mean(hs_vector(~ix_nan))) ' m elevation'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)

                hfig2 = figure('name','Surface elevation variation with potential outliers');
                plot(date_vector_J,hs_vector,'k.')
                datetick('x','mmm/yyyy')
                ylabel('Surface elevation [m]','fontsize',fontsize)
                title(['Surface elevation variation'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)

                hfig3 = figure('name','Surface temperature variation with potential outliers');
                plot(date_vector_J,Ts_vector,'k.')
                datetick('x','mmm/yyyy')
                ylabel('Surface temperature [deg]','fontsize',fontsize)
                title(['Surface Temperature variation'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)

                hfig4 = figure('name','Scaling height with potential outliers');
                plot(date_vector_J,h_scaling_vector./1000,'k.')
                datetick('x','mmm/yyyy')
                ylabel('Scaling height [km]','fontsize',fontsize)
                title(['Scaling height variation variation'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)

                hfig5 = figure('name','Mean temperature with potential outliers');
                plot(date_vector_J,T_mean_vector./1000,'k.')
                datetick('x','mmm/yyyy')
                ylabel('Mean temperature [deg]','fontsize',fontsize)
                title(['Mean temperature from elevations below ' num2str(h_max/1000) ' km'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)



                % saving of the figure when requested
                if save_fig==1
                    fprintf(['Figure with _wpo refers to with potential outliers \n'])
                    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'PI_factor_wpo.eps'];
                    set(hfig,'PaperPositionMode','auto')
                    print(hfig,'-depsc','-r150',fig_save_name)
                    clear fig_save_name

                    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'surface_elevation_wpo.eps'];
                    set(hfig2,'PaperPositionMode','auto')
                    print(hfig2,'-depsc','-r150',fig_save_name)
                    clear fig_save_name

                    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'surface_temp_wpo.eps'];
                    set(hfig3,'PaperPositionMode','auto')
                    print(hfig3,'-depsc','-r150',fig_save_name)
                    clear fig_save_name

                    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'H_scaling_factor' filesep 'H_scaling_wpo.eps'];
                    set(hfig4,'PaperPositionMode','auto')
                    print(hfig4,'-depsc','-r150',fig_save_name)
                    clear fig_save_name

                    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'Mean_temp_wpo.eps'];
                    set(hfig5,'PaperPositionMode','auto')
                    print(hfig5,'-depsc','-r150',fig_save_name)
                    clear fig_save_name

                end
                clear figure2


            end

            ix_no_outlier = [1:length(Pi_factor_vector)]';
            ix_no_outlier(ix_outlier)=[];
            % plot the same figures but remove the potential outliers
            hfig = figure('name','Pi-factor variation');
            plot(date_vector_J(ix_no_outlier),Pi_factor_vector(ix_no_outlier),'k.')
            datetick('x','mmm/yyyy')
            ylabel('PI Factor','fontsize',fontsize)
            title(['PI-factor for surface temperature at ' num2str(mean(hs_vector(~ix_nan))) ' m elevation'],'fontsize',fontsize)
            set(gca,'fontsize',fontsize)

            hfig2 = figure('name','Surface elevation variation');
            plot(date_vector_J(ix_no_outlier),hs_vector(ix_no_outlier),'k.')
            datetick('x','mmm/yyyy')
            ylabel('Surface elevation [m]','fontsize',fontsize)
            title(['Surface elevation variation'],'fontsize',fontsize)
            set(gca,'fontsize',fontsize)

            hfig3 = figure('name','Surface temperature variation');
            plot(date_vector_J(ix_no_outlier),Ts_vector(ix_no_outlier),'k.')
            datetick('x','mmm/yyyy')
            ylabel('Surface temperature [deg]','fontsize',fontsize)
            title(['Surface Temperature variation'],'fontsize',fontsize)
            set(gca,'fontsize',fontsize)

            hfig4 = figure('name','Scaling height');
            plot(date_vector_J(ix_no_outlier),h_scaling_vector(ix_no_outlier)./1000,'k.')
            datetick('x','mmm/yyyy')
            ylabel('Scaling height [km]','fontsize',fontsize)
            title(['Scaling height variation variation'],'fontsize',fontsize)
            set(gca,'fontsize',fontsize)

            hfig5 = figure('name','Mean temperature');
            plot(date_vector_J(ix_no_outlier),T_mean_vector(ix_no_outlier)./1000,'k.')
            datetick('x','mmm/yyyy')
            ylabel('Mean temperature [deg]','fontsize',fontsize)
            title(['Mean temperature from elevations below ' num2str(h_max/1000) ' km'],'fontsize',fontsize)
            set(gca,'fontsize',fontsize)


            % saving of the figure when requested
            if save_fig==1
                fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'PI_factor.eps'];
                set(hfig,'PaperPositionMode','auto')
                print(hfig,'-depsc','-r150',fig_save_name)
                clear fig_save_name

                fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'surface_elevation.eps'];
                set(hfig2,'PaperPositionMode','auto')
                print(hfig2,'-depsc','-r150',fig_save_name)
                clear fig_save_name

                fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'surface_temp.eps'];
                set(hfig3,'PaperPositionMode','auto')
                print(hfig3,'-depsc','-r150',fig_save_name)
                clear fig_save_name

                fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'H_scaling.eps'];
                set(hfig4,'PaperPositionMode','auto')
                print(hfig4,'-depsc','-r150',fig_save_name)
                clear fig_save_name

                fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'Mean_temp.eps'];
                set(hfig5,'PaperPositionMode','auto')
                print(hfig5,'-depsc','-r150',fig_save_name)
                clear fig_save_name

            end
            clear figure2



        end


        if sar_date_sounding ==1;
            % removing the outlier
            Pi_factor_vector(ix_outlier)=[];
            date_vector_J(ix_outlier)=[];
            h_scaling_vector(ix_outlier)=[];
            
            ix_ourlier = find(sum(isnan([Pi_factor_vector h_scaling_vector]),2)>=1);
            Pi_factor_vector(ix_outlier)=[];
            date_vector_J(ix_outlier)=[];
            h_scaling_vector(ix_outlier)=[];
            
            % checking those days for which the time-stamps is the SAR time
            % stamp
            ix_SAR_date = find(ifgs_dates(date_counter)==date_vector_J);
            
            if ~isempty(ix_SAR_date)
                spectrometer_scaleheight(date_counter)=nanmean(h_scaling_vector(ix_SAR_date));
                spectrometer_PIconversion(date_counter)=nanmean(Pi_factor_vector(ix_SAR_date));
            else
               if ~isempty(Pi_factor_vector) 
                    spectrometer_scaleheight(date_counter)=nanmean(h_scaling_vector);
                    spectrometer_PIconversion(date_counter)=nanmean(Pi_factor_vector);
               else
                    spectrometer_scaleheight(date_counter)=NaN;
                    spectrometer_PIconversion(date_counter)=NaN;
               end
            end
            
            if spectrometer_PIconversion(date_counter)==0
                
                keyboard
            end
        end
        
        % saving of the data
        if date_counter==size(date_start_vector,1) && sar_date_sounding ==0
            fprintf('DONE \n')
            save(save_name_final,'date_vector_J','Pi_factor_vector','hs_vector','Ts_vector','Tm_vector','h_scaling_vector','T_mean_vector','ix_outlier','h_max');
            
            spectrometer_scaleheight = nanmean(h_scaling_vector(ix_no_outlier));
            spectrometer_PIconversion = nanmean(Pi_factor_vector(ix_no_outlier));            
             
            cd(curdir)
            % updating of the parm_aps list
            setparm_aps('spectrometer_scaleheight',spectrometer_scaleheight);
            setparm_aps('spectrometer_PIconversion',spectrometer_PIconversion);
        elseif date_counter==size(date_start_vector,1) && sar_date_sounding ==1
            cd(curdir)
            
            
            % checking if some of the dates were not estimated.
            % take the mean of all existing to fill them in
            ix_nan = find(isnan(spectrometer_PIconversion)==1);
            if ~isempty(ix_nan)
                fprintf('Some dates did not have sounding data, fill the factors with the mean of the dates \n')
                datestr(ifgs_dates(ix_nan),'yyyymmdd')
                spectrometer_scaleheight(ix_nan)=nanmean(spectrometer_scaleheight);
                spectrometer_PIconversion(ix_nan) = nanmean(spectrometer_PIconversion);
            end
            
            % updating of the parm_aps list
            setparm_aps('spectrometer_scaleheight',spectrometer_scaleheight);
            setparm_aps('spectrometer_PIconversion',spectrometer_PIconversion);
            
            save(save_name_final,'spectrometer_PIconversion','spectrometer_scaleheight','ix_nan','h_max','ifgs_dates');

            
        end
        
    end

end


cd(curdir)
