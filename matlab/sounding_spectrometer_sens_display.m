function []=sounding_spectrometer_sens_display()
% function to display the powerlaw sensitivity results.
% results are save in the figures folder.
%
%     sounding_spectrometer_sens_display.m
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
% 08/2014   DB:     Save png information, update parms_aps

fontsize  =15;

% getting parameters from parms_aps file
sounding_dir = getparm_aps('sounding_dir');
error_promp_flag = getparm_aps('sounding_error_promp');
start_date = getparm_aps('sounding_start_date');
end_date = getparm_aps('sounding_end_date');
start_year_str = start_date(1:4);
end_year_str = end_date(1:4);
start_str = start_date(5:6);
end_str = end_date(5:6);
sounding_ifg_dates = getparm_aps('sounding_ifg_dates');


time_stamp = getparm_aps('sounding_time_stamp');
time_stamp_str = [];
for k=1:size(time_stamp,1)
    if k>1
        time_stamp_str = [time_stamp_str '_' time_stamp(k,:)];
    else
        time_stamp_str = [time_stamp(k,:)];
    end
end

% check in case the figure folders were deleted to remake them
if exist([sounding_dir filesep 'Spectrometer' filesep 'figures'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer' filesep 'figures']) 
end
if exist([sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor']) 
end
if exist([sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'H_scaling_factor'],'dir')~=7
    mkdir([sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor']) 
end


% loading of the data

if strcmp(sounding_ifg_dates,'y')
    % this is the PI factor for SAR dates
    load_name = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat']

%     load_name = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_SAR_dates_1month_' time_stamp_str 'Hr.mat'  ];
    load(load_name)

   
    % updating of the parm_aps list
    setparm_aps('spectrometer_scaleheight',spectrometer_scaleheight);
    setparm_aps('spectrometer_PIconversion',spectrometer_PIconversion);

    
    hfig = figure('name','Pi-factor variation','position',[ 1352         273        1073         360]);
    plot(ifgs_dates,spectrometer_PIconversion,'k.')
    if ~isempty(ix_nan)
        hold on
        plot(ifgs_dates(ix_nan),spectrometer_PIconversion(ix_nan),'ro')
        legend('Sounding data','No data, dataset average')
    end
    datetick('x','mmm/yyyy')
    ylabel('PI Factor','fontsize',fontsize)
    title(['PI-factor variation'],'fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    
    hfig2 = figure('name','Scale height variation','position',[ 1352         273        1073         360]);
    plot(ifgs_dates,spectrometer_scaleheight,'k.')
    if ~isempty(ix_nan)
        hold on
        plot(ifgs_dates(ix_nan),spectrometer_scaleheight(ix_nan),'ro')
        legend('Sounding data','No data, dataset average')
    end
    datetick('x','mmm/yyyy')
    ylabel('Scale height [km]','fontsize',fontsize)
    title(['Scale height variation'],'fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    
    
    % saving the figures
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'SAR_PI_factor.eps'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'SAR_PI_factor.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'H_scaling_factor' filesep 'SAR_H_scaling.eps'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'H_scaling_factor' filesep 'SAR_H_scaling.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name    
    
    
    
    
else
    % this is the sensitivity analysis for the PI factor
    load_name = [sounding_dir filesep 'Spectrometer' filesep 'Spectrometer_sensitivity_' time_stamp_str 'Hr_' num2str(start_year_str) start_str '_' num2str(end_year_str) end_str '.mat'];

    load(load_name,'date_vector_J','Pi_factor_vector','hs_vector','Ts_vector','Tm_vector','h_scaling_vector','T_mean_vector','ix_outlier');







    ix_no_outlier = [1:length(Pi_factor_vector)]';
    ix_no_outlier(ix_outlier)=[];
    ix_nan = isnan(hs_vector);

    % plotting the results
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
    title(['Mean temperature'],'fontsize',fontsize)
    set(gca,'fontsize',fontsize)

    % saving the figures
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'PI_factor.eps'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'PI_factor.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name

    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'surface_elevation.eps'];
    set(hfig2,'PaperPositionMode','auto')
    print(hfig2,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'PI_factor' filesep 'surface_elevation.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name

    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'surface_temp.eps'];
    set(hfig3,'PaperPositionMode','auto')
    print(hfig3,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer' filesep 'figures' filesep 'PI_factor' filesep 'surface_temp.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name

    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'H_scaling.eps'];
    set(hfig4,'PaperPositionMode','auto')
    print(hfig4,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'H_scaling.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name

    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'Mean_temp.eps'];
    set(hfig5,'PaperPositionMode','auto')
    print(hfig5,'-depsc','-r150',fig_save_name)
    fig_save_name = [sounding_dir filesep 'Spectrometer'  filesep 'figures' filesep 'H_scaling_factor' filesep 'Mean_temp.png'];
    set(hfig,'PaperPositionMode','auto')
    print(hfig,'-dpng','-r150',fig_save_name)
    clear fig_save_name
end