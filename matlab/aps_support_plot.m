function [] = aps_support_plot(technique)
% plotting the support information used for the APS correction methods
% technique flag
% when:
% 1 powerlaw
% 2 ERA-I
% 3 Meris
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
% By Bekaert David 
%

fontsize = 15;
load tca_support

%% power- law
if technique==1
    if exist('aps_p','dir')~=7
        mkdir('aps_p')
    end
    if strcmpi(getparm_aps('powerlaw_ridge_constraint'),'y')
        
        
        deminfo = powerlaw_ridges.deminfo;
        InSAR_convexhull = powerlaw_ridges.InSAR_convexhull;
        mountain_ridge = powerlaw_ridges.mountain_ridge;
        hard_ridge_flag = powerlaw_ridges.hard_ridge_flag;
        
        
        % plotting hte powerlaw ridges
        fprintf('Plotting power-law ridges \n')

        h1 = figure('name','Power-law ridges','position', [ 200   243   610   603]);
        plot( mean([deminfo.xmin deminfo.xmax]),mean([deminfo.ymax deminfo.ymin]),'g-','linewidth',2)
        hold on
        plot( mean([deminfo.xmin deminfo.xmax]),mean([deminfo.ymax deminfo.ymin]),'r-','linewidth',2)
        hold on
        plot( mean([deminfo.xmin deminfo.xmax]),mean([deminfo.ymax deminfo.ymin]),'r--','linewidth',2)

        
        imagesc([deminfo.xmin deminfo.xmax],[deminfo.ymax deminfo.ymin],deminfo.dem_smooth)   
        cc = colorbar;
        view(0,90)
        axis equal
        axis tight   
        colormap(flipud(gray))
        axis xy
        box on

        hold on
        plot(InSAR_convexhull(:,1),InSAR_convexhull(:,2),'g-','linewidth',2)

        % plotting the line on the figure
        for counter=1:length(mountain_ridge)
            if hard_ridge_flag(counter)==1
                plot(mountain_ridge{counter}(:,1),mountain_ridge{counter}(:,2),'r-','linewidth',2)
            else
                plot(mountain_ridge{counter}(:,1),mountain_ridge{counter}(:,2),'r--','linewidth',2)
            end
        end
        box on
        xlabel(cc,'[m]','fontsize',fontsize)
        ylabel(cc,'Topography','fontsize',fontsize)
        legend('Location','NorthOutside','InSAR region','Hard ridge','Subregion boundary')
        title('Ridge definition','fontsize',fontsize)
        set(gca,'fontsize',fontsize)
        print(h1,'-dpng',['aps_p/mountain_ridges.png'])
        print(h1,'-depsc',['aps_p/mountain_ridges.eps'])

    else
       % plotting the power-law windows 
        fprintf('Plotting power-law windows \n')
        
        % see if the real dem is available otherzie use the poitn heights
        dem_flag=1;
        if exist('powerlaw_ridges','var')==1 
            deminfo = powerlaw_ridges.deminfo;
        elseif exist('era','var')==1 
            deminfo = era.deminfo;
        else
            hgt_file = getparm_aps('hgt_matfile');
            ll_matfile = getparm_aps('ll_matfile');
    
            hgt =  load(hgt_file);
            ll = load(ll_matfile);
            if strcmpi(getparm_aps('stamps_processed'),'y')
                hgt = hgt.hgt;
                ll = ll.lonlat;
                dem_flag=0;
            end
        end
        % making sure the dem and the other data have the same origin in
        % longitude
        if dem_flag== 1
            if deminfo.xmax-max(InSAR_convexhull(:,1))>100
                deminfo.xmax = deminfo.xmax-360;
                deminfo.xmin = deminfo.xmin-360;
            end
        end
        
        % laoding of the data
        window_box_ll= powerlaw_windows.window_box_ll;
        window_box_center_ll = powerlaw_windows.window_box_center_ll;
        
        
        % plotting the results
        h1 = figure('name','Power-law windows','position', [ 200   243   610   603]);
        % for the legend
        plot(window_box_ll{round(0.5*(length(window_box_ll)))}(:,1),window_box_ll{round(0.5*(length(window_box_ll)))}(:,2),'r-','linewidth',2)
        hold on
        plot(window_box_center_ll(1,1),window_box_center_ll(1,2),'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
        hold on
        plot(InSAR_convexhull(:,1),InSAR_convexhull(:,2),'g-','linewidth',2)
        hold on

        if dem_flag==1
             % plot the topogrpahy
            imagesc([deminfo.xmin deminfo.xmax],[deminfo.ymax deminfo.ymin],deminfo.dem)  
            view(0,90)
            axis xy
        else
            scatter3(ll(:,1),ll(:,2),hgt,3,hgt,'filled');
            view(0,90)
        end
        cc=colorbar;
        colormap(flipud(gray))
        xlabel(cc,'[m]','fontsize',fontsize)
        ylabel(cc,'Topography','fontsize',fontsize)
        axis equal
        axis tight 
        
        hold on
        plot3(window_box_ll{round(0.5*(length(window_box_ll)))}(:,1),window_box_ll{round(0.5*(length(window_box_ll)))}(:,2),99999.*ones(size(window_box_ll{1}(:,2))),'r-','linewidth',2)
        hold on
        plot3(window_box_center_ll(:,1),window_box_center_ll(:,2),10000000.*ones(size(window_box_center_ll(:,2))),'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
        hold on
        plot3(InSAR_convexhull(:,1),InSAR_convexhull(:,2),99999.*ones(size(InSAR_convexhull(:,2))),'g-','linewidth',2)
       
  
    
        legend('Location','NorthOutside','Window example','Window centers','InSAR region')
        title('Window definition','fontsize',fontsize)
        set(gca,'fontsize',fontsize)

        grid off
        box on
        
        print(h1,'-dpng',['aps_p/window_locations.png'])
        print(h1,'-depsc',['aps_p/window_locations.eps'])
        
    end
 
%% ERA interim
elseif technique==2
    if exist('aps_e','dir')~=7
        mkdir('aps_e')
    end
    
    % loading the ERA data
    deminfo = era.deminfo;
    era_lonlat = era.era_lonlat;
    
    % check if the InSAR outline exist, if so plot it too
    if exist('InSAR_convexhull','var')==1;
        insar_flag =1;
    else
        insar_flag=0;
    end
    if insar_flag== 1
        if deminfo.xmax-max(InSAR_convexhull(:,1))>100
            InSAR_convexhull(:,1) = InSAR_convexhull(:,1)+360;
        end
    end
    
    % plotting the figure
    hfig = figure('name','ERA points and InSAR region','position', [ 200   243   610   603]);
    % for the legend generation
    plot(mean([deminfo.xmin deminfo.xmax]),mean([deminfo.ymax deminfo.ymin]),'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
    hold on
    % check if the insar convex hull can be plotted
    if insar_flag==1
        plot(mean([deminfo.xmin deminfo.xmax]),mean([deminfo.ymax deminfo.ymin]),'g-','linewidth',2)
    end
    % plot the topogrpahy
    imagesc([deminfo.xmin deminfo.xmax],[deminfo.ymax deminfo.ymin],deminfo.dem)  
    cc=colorbar;
    view(0,90)
    colormap(flipud(gray))
    axis xy
    xlabel(cc,'[m]','fontsize',fontsize)
    ylabel(cc,'Topography','fontsize',fontsize)
    hold on
    plot(era_lonlat(:,1),era_lonlat(:,2),'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
    hold on
    if insar_flag==1
        plot(InSAR_convexhull(:,1),InSAR_convexhull(:,2),'g-','linewidth',2)
        legend('location','northoutside','Used ERA-I locations','InSAR region')
    else
        legend('location','northoutside','Used ERA-I locations')
    end
    title('ERA points distribution','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    axis equal
    axis tight 
    
    print(hfig,'-dpng',['aps_e/era_datapoints.png'])
    print(hfig,'-depsc',['aps_e/era_datapoints.eps'])

    
elseif technique==3 % MERIS
    load psver
    ps = load(['ps' num2str(psver) '.mat']);
    stamps_processed = getparm_aps('stamps_processed');
    small_baseline_flag  = getparm('small_baseline_flag');
    if strcmpi(small_baseline_flag,'y')
        load('tca_sb2.mat','ph_tropo_meris');
        % for baselines you need ps of PS directory 
        bperp = load(['..' filesep 'ps1.mat'],'bperp');
        bperp = bperp.bperp;
        bperp = bperp(ps.ifgday_ix);
        dates=ps.ifgday;   
        ix_ifgs_keep = 1:ps.n_ifg;
    else
        load('tca2.mat','ph_tropo_meris');
        dates = [ps.day   repmat(ps.master_day,length(ps.day),1)]; 
        bperp = [ps.bperp zeros([length(ps.day) 1])];
        ix_ifgs_keep = 1:ps.n_ifg;
    end

    
    % plot the SB network and show for whcih connection there is meris data
    if strcmp(stamps_processed,'y')

        % remove that that have not been unwrapped before
        ix_dropped = getparm('drop_ifg');
        
        % interferograms with a meris correction
        temp = ix_ifgs_keep;
        temp(ix_dropped)=[];
        ifgs_good = find(sum(ph_tropo_meris(:,temp)~=0,1)~=0)
        
        
        if strcmp(small_baseline_flag,'n')
            ix_dropped = unique([ix_dropped ps.master_ix]);
        end
        ix_ifgs_keep(ix_dropped)=[];
        bperp(ix_dropped,:)=[];
        dates(ix_dropped,:)=[];

        % keep the original network for reference
        dates_all = dates;
        bperp_all = bperp;



        % plotting the network of the interferograms used in the RMSE computation
        h_baselineplot = figure('name','Processed network');
        % plotting all ifgs that were considered
        for ifgs_counter=1:size(bperp_all,1)
            hold on
            plot([dates_all(ifgs_counter,1)  dates_all(ifgs_counter,2)],  [bperp_all(ifgs_counter,1) bperp_all(ifgs_counter,2)] ,'k-','linewidth',1) 
        end
        


        % plotting the network for which we have an APS correction 
        for ifgs_counter=1:length(ifgs_good)
            ifgs_counter
            hold on
            plot([dates(ifgs_good(ifgs_counter),1)  dates(ifgs_good(ifgs_counter),2)],  [bperp(ifgs_good(ifgs_counter),1) bperp(ifgs_good(ifgs_counter),2)] ,'k-','linewidth',2) 
            text(mean(dates(ifgs_good(ifgs_counter),:)),mean(bperp(ifgs_good(ifgs_counter),:))+20,num2str(ifgs_good(ifgs_counter)),'fontsize',fontsize-2,'backgroundColor',[1 1 1],'Margin',0.01)
        end
        hold on
        % [dates_unique,ix_unique] = unique(dates(ifgs_good,:));
        % bperp_unique = bperp(ifgs_good,:);
        % bperp_unique = bperp_unique(ix_unique);
        [dates_unique,ix_unique] = unique(dates_all);
        bperp_unique = bperp_all;
        bperp_unique = bperp_unique(ix_unique);
        plot(dates_unique,bperp_unique,'ko','markerfacecolor','r','markersize',7)
        clear dates_unique ix_unique bperp_unique
        hold on
        % set(gca,'XTick',dates_num)
        datetick('x','mmm yy')
        set(gca,'fontsize',fontsize)
        box on
        ylabel('Bperp [m]','fontsize',fontsize)
    end

    
    
end