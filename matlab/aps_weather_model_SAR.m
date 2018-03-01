function [] = aps_weather_model_SAR(model_type)
% [] = aps_weather_model_SAR()
% Script to load weather model data and compute SAR delays
% The DEM file can be .grd file or other. If the latter, the DEM should have an asociated
% ".rsc" file, with the same filename as the DEM. The ".rsc" files should
% contain a WIDTH, LENGTH, X_FIRST, Y_FIRST, X_STEP, Y_STEP and optional a 
% FORMAT string. The weather model data is assumed to be structured in date_before(d,:) folders. 
%
%
% INPUTS:
% demfile               Full path to the DEM file. The DEM needs to me in meters.
% xlims                 Limits in the x-direction, either in degrees
% ylims                 Limits in the y-direction, either in degrees
% demnull               The value for no DEM data, default is -32768.
% smpres                The output resolution, either in degrees
%                       Units needs to be consistend with xlims and ylims.
%
% OUTPUTS:
% It will give the computed ZENITH hydrostatic and wet delay map in cm for the selected region. 
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
%
% Modified from Richard Walters - Oxford / Leeds University 2012-2013
% Modifcations:
% DB 	10/2013		Convert script to a function, integrate in aps_toolbox, 
%                   add syntax. Retrieve parameters automatic.
% DB    10/2013     Output the combined delay map.
% DB    10/2013     Increased processing efficiency
% DB    10/2013     add compatibility with Doris construct_dem.sh script
% DB    10/2013     Bug fix incase the the netcdf is split over two days
% DB    11/2013     Change the computation f saturated pressure to be
%                   identical witt how it was computed in ERA-I model.
% DB    11/2013     Update filenaming to Hydrostatic as its not dry delay.
%                   Compute the hydrostatic delay directly from int(P/T)
% DB    01/2014     Fix bug for the ERA-I longitude boundary, being set to small. 
% DB    01/2014     Allow for tiled ERA-I data
% HB-DB 02/2014     Fix to include ECMWF with scaling and offset
% DB    02/2014     Include a parameter for ECMWF data website
% DB    03/2014     Include an extra check for the DEM grd file and the
%                   selected crop
% DB 	03/2014     Save individual time stamp SAR date_before(d,:) processing
% DB    07/2014     Redefine meris_lat(lon)_range to region_lat(lon)_range
% DB    08/2014     Check if netcdf exist (ERA-I is not an operation service)
% DB    11/2014     Include option to output the 3D delays for each date_before(d,:)
% HW    02/2015     check dem file format automatically
% DB 	06/2015     Bug fix for plotting the support information 
% DB    06/2015     Fix typo in error message
% DB    11/2015     Branch of the DEM in to seperate function same for all codes
% DB    11/2015     Add multicore option from matlab
% DB 	02/2016     Close netcdf files
% DB    04/2016     Branch into weather model script and include merra too
% SSS   04/2016     Clear variables such memory need is reduced
% DB    07/2016     redefine hydrostatic delay to be based on surface pressure.
% DB    08/2016     Uncomment a keyboard
% DB    07/2016     expand to include ERA5 test model data
% KM    02/2018     expand to include NARR model data

fig_test = 0;           % when 1 show the dem as debug figure
save_3D_delays = 0;     % When 1 saves the tropopsheric delays for each x,y and with height

if nargin<1
    error('Give at least the model_type: era, era5, narr, merra, or merra2')
end
% change to lower caps for saving and filename generation consistency
model_type = lower(model_type);

%% Constants
% Parameters from Jolivet et al 2011, GRL
Rd = 287.05;            % [j/kg/K] specific gas constants dry air
Rv = 461.495;           % [j/kg/K] water vapour
k1 = 0.776;             % [K/Pa]
k2 = 0.716;             % [K/Pa]
k3 = 3.75e3;            % [K^2/Pa]
coeff.Rd = Rd;
coeff.Rv = Rv;

%% Defaults
zref = 15000;       % zref for integration- tropopause
zincr = 15;         % z increment spacing for vertical interpolation before integration
vertres = 5;        % vertical resolution for comparing dem to vert profiles of delay

%% output file names
% output file for the DEM and look angles
smpdem = 'dem_smp.xyz';

% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed',1);

%%% defaults for the weather models. If not applicable it will be changed below for the specific model.
if strcmp(model_type,'narr')
    timelist_model = ['0000' ;'0300'; '0600' ; '0900'; '1200' ;'1500'; '1800' ;'2100'; '0000'];
    model_lag = 8*7;    % days
else
    timelist_model= ['0000' ; '0600' ; '1200' ; '1800' ; '0000'];       % the time interval the model is outputed
    model_lag = 0; % ? check lags
end
era_data_type = [];                                                 % the weather model data type for ERA only.


%%% Updating specific weather model information
if strcmpi(model_type,'era')
    weather_model_datapath = getparm_aps('era_datapath',1);
    era_data_type = getparm_aps('era_data_type');         % the datatype of the model either BADC or ERA
elseif strcmpi(model_type,'merra') || strcmpi(model_type,'merra2')
    weather_model_datapath = getparm_aps('merra_datapath',1); 
elseif strcmpi(model_type,'narr') 
    weather_model_datapath = getparm_aps('narr_datapath',1); 

else
    error(['weather model type not supported, either: wrf, era, narr, merra for now'])
end

lambda = getparm_aps('lambda',1)*100;                       % radar wavelength in cm
demfile = getparm_aps('demfile',1);
demnull = getparm_aps('dem_null',1);
UTC_sat =  getparm_aps('UTC_sat',1);
ifgday_matfile = getparm_aps('ifgday_matfile',1);
ifgs_dates = load(ifgday_matfile);


% loading the data
if strcmp(stamps_processed,'y')
   dates = ifgs_dates.day;
   load psver
   fprintf('Stamps processed structure \n')
else
    psver = 2;
    ifgs_dates = ifgs_dates.ifgday;
    dates = reshape(ifgs_dates,[],1);
    dates = unique(dates);
    dates = datenum(num2str(dates),'yyyymmdd');
end

% getting the dates
n_dates = length(dates);

%% Compute and resample DEM
[dem,xmin,xmax,ymin,ymax,smpres,nncols,nnrows] = get_DEM;

% the region which is cropped from the ERA data and used to make the interpolation.
% Should be  larger than the region to which the delay is computed
lonmin = floor(xmin)-1;
lonmax= ceil(xmax)+1;
latmin = floor(ymin)-1;
latmax = ceil(ymax)+1;

% setting the maximum height of the DEM to limit the range at which ERA-I
% needs to be interpolated to
maxdem = ceil(max(max(dem))/100)*100+50; % max height of dem in metres

fprintf(['Interpolate to a maximum dem height of ', num2str(maxdem) ,' m\n'])


%% Compute based on Satellite pass which weather model outputs that will be used
[time_before,time_after, date_before, date_after,f_before,f_after] = aps_weather_model_times(timelist_model,dates,UTC_sat,model_lag);

%% generating a file 
[modelfile_before,modelfile_after] = aps_weather_model_filenames(model_type,time_before,time_after,date_before, date_after,weather_model_datapath);


%% performing the calucluation for each date 
for d = 1:n_dates
 
    % the save filenames
    outfile = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZWD.xyz'];
    hydroutfile = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZHD.xyz']; 
     
    no_data = 0;
    for kk = 1:2
        if kk == 1
            file = modelfile_before(d,:);
        end
        if kk == 2
            file = modelfile_after(d,:);
        end
        
        
        % check
        if exist(file,'file')~=2
            no_data = no_data+1;
        else
               
            % loading the weather model data
            if  strcmpi(model_type,'era')
                 [ Temp,e,Geopot,P,longrid,latgrid,xx,yy,lon0360_flag] =  aps_load_era(file,era_data_type) ;
            elseif  strcmpi(model_type,'era5')
                 [ Temp,e,Geopot,P,longrid,latgrid,xx,yy,lon0360_flag] =  aps_load_era(file,era_data_type) ;
            elseif strcmpi(model_type,'merra') || strcmpi(model_type,'merra2')
                 [ Temp,e,Geopot,P,longrid,latgrid,xx,yy,lon0360_flag] =  aps_load_merra(file,model_type,coeff) ;
            elseif strcmpi(model_type,'narr')
                 [ Temp,e,Geopot,P,longrid,latgrid,xx,yy,lon0360_flag] =  aps_load_narr(file,model_type,coeff) ;
            end
%             longrid=longrid';latgrid=latgrid';
            % verify and cope with NAN's
            [Temp,e,Geopot,P,longrid,latgrid] =  aps_weather_model_nan_check(Temp,e,Geopot,P,longrid,latgrid) ;

            % define weather model grid nodes
            latlist = reshape(latgrid(:,:,1),[],1);
            lonlist = reshape(longrid(:,:,1),[],1);
            xlist = reshape(xx,[],1);
            ylist = reshape(yy,[],1);

            % Limit weather model to those grid points around the user-defined InSAR box
            % Getting the weather model resolution
            lat_res = abs(diff(unique(latgrid)))*1.5;
            lat_res = lat_res(1);
            lon_res = abs(diff(unique(longrid)))*1.5;
            lon_res = lon_res(1);

            % make sure the weather model and the lonlat specified by user
            % is defined in the same way. Weatehr models can be [0 360]
            % degrees. User can be [-180 180] unlikely [0 360]
            if strcmpi(lon0360_flag,'y')
                % fprintf('This needs to be defined better and some checks are needed \n')
                % not all pixels should be shifted.
  
                if xmin<0
                   xmin= xmin+360; 
                end
                if xmax<0
                   xmax = xmax +360; 
                end
            end

            % generation a plot
            if fig_test ==1 & d==1 & kk==1
                fontsize = 15;
                hfig = figure('name',[model_type ' weather model nodes and InSAR region'],'position', [ 200   243   610   603]);
                % for the legend generation
                plot(mean([xmin xmax]),mean([ymax ymin]),'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
                hold on
                plot(mean([xmin xmax]),mean([ymax ymin]),'r.')
                hold on
                plot(mean([xmin xmax]),mean([ymax ymin]),'r-','linewidth',2)

                imagesc([xmin xmax],[ymax ymin],dem)  
                cc=colorbar;
                view(0,90)

                axis xy
                ylabel(cc,'Topography [m]','fontsize',fontsize)
                hold on
                plot(lonlist,latlist,'wo','markeredgecolor','k','markerfacecolor','w','markersize',15)
                hold on
                plot([xmin xmin xmax xmax xmin],[ymin ymax ymax ymin ymin],'r-','linewidth',2)

                title([model_type ' weather model nodes'],'fontsize',fontsize)
                set(gca,'fontsize',fontsize)
                axis equal
                axis tight 

                legend('location','northoutside',[model_type ' weather model nodes'],['Used ' model_type ' weather model nodes'],'InSAR box')
            end

            % making the weather model grid slightly larger than the InSAR
            % bounding box. Will use the weather model resolution for this
            % to make sure an extra grid node is included.
            ix = find(ymin-lat_res<= latlist & latlist <= ymax+lat_res  & xmin-lon_res<= lonlist & lonlist<= xmax+lon_res) ;
            xlist = xlist(ix);
            ylist = ylist(ix);
            latlist = latlist(ix);
            lonlist = lonlist(ix);
            numy = length(unique(latlist));
            numx = length(unique(lonlist));
            ulatlist = unique(latlist);
            ulonlist = unique(lonlist);
            uxlist = unique(xlist);
            uylist = unique(ylist);

            if fig_test ==1 & d==1 & kk==1
                figure(hfig)
                hold on
                plot(lonlist,latlist,'r.')
                xlim([xmin-4*lon_res  xmax+4*lon_res])
                ylim([ymin-4*lat_res  ymax+4*lat_res])
                str='';
                while ~strcmpi(str,'y') && ~strcmpi(str,'n') 
                    fprintf(['Do the red nodes extend (just) outside the red InSAR rectangle? \n'])  
                    str = input('Continue? [y: for yes, n: no] \n','s');
                end
                if strcmpi(str,'n')
                    error('Check your lon lat crop, otherwize extend area of downlaoded weather model data.')
                end
            end
            
            clear lon_res lat_res ylist xlist %%%SSS 4/16

            % saving the information for support plotting 
            eval([model_type '.' model_type '_lonlat =[lonlist latlist];']);                            % era.era_lonlat = [lonlist latlist];
            eval([model_type '.region =[[xmin xmin xmax xmax xmin]'' [ymin ymax ymax ymin ymin]''];']); % era.region = [[xmin xmin xmax xmax xmin]' [ymin ymax ymax ymin ymin]'];
            deminfo.xmin = xmin;
            deminfo.xmax = xmax;
            deminfo.ymax = ymax;
            deminfo.ymin = ymin;
            deminfo.dem = dem;
            eval([model_type '.deminfo =deminfo;']);                                                    % era.deminfo = deminfo;
            clear deminfo
            % checking if the file already exist. Yes append otherwiuze create it
            if exist('tca_support.mat','file')==2
                eval(['save(''tca_support.mat'',''-append'',''' model_type ''');']);                    % save('tca_support.mat','-append','era')
            else
                eval(['save(''tca_support.mat'',''' model_type ''');']);                                % save('tca_support.mat','era')        
            end
%             if fig_test ==1 & d==1 & kk==1
%                     if exist(['aps_' model_type],'dir')~=7
%                         mkdir(['aps_' model_type]);
%                     end
%                     print(hfig,'-dpng',['aps_' model_type  filesep model_type '_datapoints.png'])
%                     print(hfig,'-depsc',['aps_' model_type filesep model_type '_datapoints.eps'])
%             end
            eval(['clear ' model_type]);
         

            % Convert Geopotential to Geopotential Height and then to Geometric Height
            g0 = 9.80665;
            % Calculate Geopotential Height, H
            H = Geopot./g0;

            % map of g with latitude
            g = 9.80616.*(1 - 0.002637.*cosd(2.*latgrid) + 0.0000059.*(cosd(2.*latgrid)).^2);
            % map of Re with latitude
            Rmax = 6378137; 
            Rmin = 6356752;
            Re = sqrt(1./(((cosd(latgrid).^2)./Rmax^2) + ((sind(latgrid).^2)./Rmin^2)));

            % Calculate Geometric Height, Z
            Z = (H.*Re)./(g/g0.*Re - H);
            
            % Find middle of scene to work out glocal and Rlocal for use later
            midx = round(mean(uxlist));
            midy = round(mean(uylist));
            glocal = g(midy,midx,1);
            Rlocal = Re(midy,midx,1);

            cdslices = maxdem/vertres +1;
            cdstack = zeros(numy,numx,cdslices);
            cdstack_dry = zeros(numy,numx,cdslices);
            cdstack_wet = zeros(numy,numx,cdslices);

            XI=(0:zincr:zref)';
            gh = glocal.*(Rlocal./(Rlocal+XI)).^2; %gravity with height for XI height range

            % Interpolate Temp P and e from 0:20:15000 m
            % then integrate using trapz to estimate delay as function of height
            for i = 1:numx;
                for j = 1:numy;
                    xn = uxlist(i);
                    yn = uylist(j);

                    %interpolate at res zincr before doing integration
                    X=double(squeeze(Z(yn,xn,:)));
                    Ye=double(squeeze(e(yn,xn,:)));
                    YeI = interp1(X,Ye,XI,'spline')*100; %change from hPa to Pa

                    YP=double(squeeze(P(yn,xn,:)));
                    YPI = interp1(X,YP,XI,'spline')*100; %change from hPa to Pa

                    YT=double(squeeze(Temp(yn,xn,:)));
                    YTI = interp1(X,YT,XI,'spline');

                    tmp1 = ((k2-(Rd*k1/Rv)).*YeI./YTI + k3.*YeI./(YTI.^2));
                    Lw = (10^-6).*-1*flipud(cumtrapz(flipud(XI),flipud(tmp1)));
                    % L is zenith delay one way in metres
                    % tmp2 = k1.*YPI./YTI;
                    %Ld = (10^-6).*-1*flipud(cumtrapz(flipud(XI),flipud(tmp2)));             % This is using P/T expression (Hanssen, 2001)
                    gm = glocal; 
                    Ld = (10^-6).*((k1*Rd/gm).*(YPI - YPI(zref/zincr +1)));                 % This is P0 expression (Hanssen, 2001)

                    % Interpolate important part (i.e. total delay at elevations
                    % less than maxdem) at high res i.e. vertres, and put in cdstack.
                    cdI=(0:vertres:maxdem)';
                    LdI=interp1(XI,Ld,cdI,'spline');
                    LwI=interp1(XI,Lw,cdI,'spline');

                    %cdstack(j,i,:) = LsI;
                    cdstack_dry(j,i,:) = LdI;
                    cdstack_wet(j,i,:) = LwI;
                    
                    if save_3D_delays==1               
                        cdI3D=(0:100:max(XI))';
                        LdI3D=interp1(XI,Ld,cdI3D,'spline');
                        LwI3D=interp1(XI,Lw,cdI3D,'spline');

                        %cdstack(j,i,:) = LsI;
                        cdstack_dry3D(j,i,:) = LdI3D;
                        cdstack_wet3D(j,i,:) = LwI3D;
                        clear LdI3D LwI3D
                    end
                end
            end
            clear uxlist uylist ulonlist ulatlist %%%SSS 4/16            

            % Interpolate each cdstack layer onto a grid given by the DEM extents
            % in UTM m.
            xsmpres = (xmax-xmin)/nncols;
            ysmpres = (ymax-ymin)/nnrows;
            [xi,yi] = meshgrid(xmin+0.5*xsmpres:xsmpres:xmax-0.5*xsmpres,ymin+0.5*ysmpres:ysmpres:ymax-0.5*ysmpres);

            ix_temp = diff(lonlist);
            ix_temp = find(ix_temp~=0);
            ix_temp = ix_temp(1);
            lonlist_matrix = reshape(lonlist,ix_temp,[]) ;
            latlist_matrix = reshape(latlist,ix_temp,[]) ;
            clear lonlist latlist %%%SSS 4/16
            
            % saving the outputs
            if save_3D_delays==1               
                if kk==1
                    clear hgt lon lat dry1 dry2 wet1 dem_temp hgt_topo %%%SSS 4/16
                    hgt = cdI3D;
                    lon = [lonlist_matrix];
                    lat = latlist_matrix;
                    dry1 = cdstack_dry3D;
                    wet1 = cdstack_wet3D;
                                              
                    % also give the station topography
                    dem_temp = dem;
                    dem_temp(isnan(dem_temp))=0;
                    hgt_topo = griddata(xi,yi,dem_temp,lon,lat,'linear'); 
                    clear dem_temp
                    
                    
                elseif kk==2
                    clear dry2 wet2
                    dry2 = cdstack_dry3D;
                    wet2 = cdstack_wet3D;
                end
            end
            
            clear cdstack_interp_dry %%%SSS 4/16
            cdstack_interp_dry = zeros(nnrows,nncols,cdslices);
            parfor n = 1:cdslices
                newslice = interp2(lonlist_matrix,latlist_matrix,cdstack_dry(:,:,n),xi,yi,'linear');          
                cdstack_interp_dry(:,:,n)= flipud(newslice); % flipud due to ypositive downpage for matlab figs, and ypositive uppage for utmy
            end     

            clear cdstack_interp_wet %%%SSS 4/16 
            cdstack_interp_wet = zeros(nnrows,nncols,cdslices);
            parfor n = 1:cdslices 
               newslice = interp2(lonlist_matrix,latlist_matrix,cdstack_wet(:,:,n),xi,yi,'linear');   
               cdstack_interp_wet(:,:,n)= flipud(newslice); % flipud due to ypositive downpage for matlab figs, and ypositive uppage for utmy
            end
            clear lonlist_matrix latlist_matrix %%%SSS 4/16
            % keeping the coordinates in the same grid as the data
            xi = flipud(xi);
            yi = flipud(yi);


            % Pull out delays from cdstack layers that match dem heights
            clear wetcorrection hydrcorrection rounddem %%%SSS 4/16
            wetcorrection = ones(nnrows,nncols);
            hydrcorrection = ones(nnrows,nncols);
            rounddem = round(dem/vertres);
            rounddem(dem < 0)=0;

            % cope with the case that NaN are present. This can be sea-level
            rounddem(isnan(dem))=0;

            for i=1:nnrows
                for j=1:nncols
                    wetcorrection(i,j) = cdstack_interp_wet(i,j,rounddem(i,j)+1);
                end
            end

            for i=1:nnrows
                for j=1:nncols
                    hydrcorrection(i,j) = cdstack_interp_dry(i,j,rounddem(i,j)+1);
                end
            end

            if kk==1
                wetcorrection1 = wetcorrection;
                drycorrection1 = hydrcorrection;
            end
            if kk==2
                wetcorrection2 = wetcorrection;
                drycorrection2 = hydrcorrection;
            end
            clear wetcorrection hydrcorrection

        end
    end

    if sum(no_data)==0
        % note that this is a one way Zenith delay and not a slant delay. 
        % Units are in cm

        % saving individual estimates based on the time-stamp
        outfile_wet_before = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZWD_before.xyz'];
        outfile_wet_after = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZWD_after.xyz'];
        outfile_dry_before = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZHD_before.xyz'];
        outfile_dry_after = [weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_ZHD_after.xyz'];


        fidwet_before = fopen(outfile_wet_before,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(wetcorrection1,[],1)]';
        tally = fwrite(fidwet_before,data_write,'double');
        fclose(fidwet_before);
        clear data_write tally %%%SSS 4/16
        fidwet_after = fopen(outfile_wet_after,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(wetcorrection2,[],1)]';
        tally = fwrite(fidwet_after,data_write,'double');
        fclose(fidwet_after);
        clear data_write tally %%%SSS 4/16
        fiddry_before = fopen(outfile_dry_before,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(drycorrection1,[],1)]';
        tally = fwrite(fiddry_before,data_write,'double');
        fclose(fiddry_before); 
        clear data_write tally %%%SSS 4/16
        fiddry_after = fopen(outfile_dry_after,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(drycorrection2,[],1)]';
        tally = fwrite(fiddry_after,data_write,'double');
        fclose(fiddry_after);
        clear data_write tally %%%SSS 4/16

        % saving the outputs
        if save_3D_delays==1   
            wet = wet1*f_before(d) +  wet2*f_after;
            dry = dry1*f_before(d) +  dry2*f_after;
            save([weather_model_datapath filesep date_before(d,:) filesep date_before(d,:) '_3D.mat'],'dry','wet','hgt','lon','lat','hgt_topo')
            clear wet dry hgt dry1 dry2 wet1 wet2 %%%SSS 4/16
        end
                
        % Output wet correction
        wetcorrection = wetcorrection1*f_before(d) +  wetcorrection2*f_after(d);
        clear wetcorrection1 wetcorrection2
        wetcorrection = wetcorrection*100;                 % delay in [cm]
        fid = fopen(outfile,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(wetcorrection,[],1)]';
        tally = fwrite(fid,data_write,'double');
        fclose(fid);
        clear data_write tally wetcorrection %%%SSS 4/16

        % Output hydrostatic correction
        hydrcorrection = drycorrection1*f_before(d) +  drycorrection2*f_after(d);
        clear drycorrection1 drycorrection2
        hydrcorrection = hydrcorrection*100;                 % delay in [cm]
        fid = fopen(hydroutfile,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(hydrcorrection,[],1)]';
        tally = fwrite(fid,data_write,'double');
        fclose(fid);
        clear data_write tally hydrcorrection %%%SSS 4/16

        fprintf([num2str(d) ' completed out of ' num2str(n_dates) '\n' ])
        
    else
        fprintf([num2str(d) ' completed out of ' num2str(n_dates) ' (NO DATA)\n' ])
    end
end


