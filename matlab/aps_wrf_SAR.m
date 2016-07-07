function [] = aps_wrf_SAR()
% [] = aps_wrf_SAR()
% Script to load wrf processed data, on 37 pressure levels. Load data for area of interest, 
% and calculate delays maps for each specified date.
% The DEM file inputed should have an asociated
% ".rsc" file, with the same filename as the DEM. The ".rsc" files should
% contain a WIDTH, LENGTH, X_FIRST, Y_FIRST, X_STEP, Y_STEP and optional a 
% FORMAT string. The meris data is assumed to be structured in date folders. 
% The batchfile contains the full path to the meris files in these folders. 
% Note that the first line of the batchfile should read "files".
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
% It will give the computed ZENITH dry and wet delay map in cm for the selected region. 
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
% Richard Walters/ Bekaert David - Oxford / Leeds University 2012-2013
% modifications:
% 19/11/2013    DB      Change the computation of the saturated vapor
%                       pressure over ice and water
% 19/11/2013    DB      Change the computation of the hydrostatic delay to
%                       be directly computed from pressure and temperature
% 20/11/2013    DB      Check if wrf data exist to avoid crashing 
% 21/11/2013    DB      Interpolate H,Temp, Hum spatially to cope with
%                       coarser land mask compared to the DEM
% 02/01/2014	DB      fix for NaN value over sea
% 04/08/2014    DB      Redefine meris_lat(lon)_range to region_lat(lon)_range
% 24/02/2015    DB      Save support information for WRF
% 27/09/2015    DB      Separate the DEM loading, ask for number of domains.
% 28/09/2015    DB      Also save the delay computation for a grid or z to include Pirate compartibility.
% 28/09/2015    DB      Include multi-core from matlab
% 07/07/2016    DB      Redefine hydrostatic delay to be based on surface pressure.
 
save_complete=0;        % save support information when 0


%% Constants
% Parameters from Jolivet et al 2011, GRL
Rd = 287.05;            % [j/kg/K] specific gas constants dry air
Rv = 461.495;           % [j/kg/K] water vapour
k1 = 0.776;             % [K/Pa]
k2 = 0.716;             % [K/Pa]
k3 = 3.75e3;            % [K^2/Pa]


%% Defaults
zref = 15000;       % zref for integration- tropopause
zincr = 15;         % z increment spacing for vertical interpolation before integration
vertres = 5;        % vertical resolution for comparing dem to vert profiles of delay


%% get the number of domains
repeat=1;
while repeat==1
    action_flag= str2num(input('How many domains do you have including the parent? ','s'));
    if isnumeric(action_flag)
        n_domains = action_flag;
        repeat=0;
    end
end


%% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed',1);
ll_matfile = getparm_aps('ll_matfile',1);
wrf_datapath = getparm_aps('wrf_datapath',1);
datestructure = 'yyyymmdd';                               % assumed date structure for ERA-I
UTC_sat =  getparm_aps('UTC_sat',1);

% loading the data
if strcmp(stamps_processed,'y')
   ps = load(ll_matfile);
   dates = ps.day;
   load psver
   fprintf('Stamps processed structure \n')
else
    psver = 2;
    ifgday_matfile = getparm_aps('ifgday_matfile',1);
    ifgs_dates = load(ifgday_matfile);
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


% setting the maximum height of the DEM to limit the range at whcih ERA-!
% needs to be interpolated to
maxdem = ceil(max(max(dem))/100)*100+50; % max height of dem in metres

fprintf(['Interpolate to a maximum dem height of ', num2str(maxdem) ,' m\n'])




%% performing the calucluation for each  date

for d = 1:n_dates
    date = datestr(dates(d),datestructure);
    date_str_file = datestr(dates(d),'yyyy-mm-dd');
    
    % the netcdf files
    infile = [wrf_datapath filesep date filesep 'wrfplev_d0' num2str(n_domains) '_' date_str_file '_' UTC_sat ':00'];
    infile2 = [wrf_datapath filesep date filesep 'wrfout_d0' num2str(n_domains) '_' date_str_file '_' UTC_sat ':00'];

    if exist(infile,'file')==2
        % the output file names
        outfile = [wrf_datapath filesep date filesep date '_ZWD.xyz'];
        hydroutfile = [wrf_datapath filesep date filesep date '_ZHD.xyz']; 
        outfile2 = [wrf_datapath filesep date filesep date '_ZWD.z'];
        hydroutfile2 = [wrf_datapath filesep date filesep date '_ZHD.z'];         
        
        % getting the information from the output file
        ncid = netcdf.open(infile,'NC_NOWRITE');
        [numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);
        [dimname, dimlen] = netcdf.inqDim(ncid,0); % dimensions
        % 3 = pressure*37 
        % 6 = temperature 

        % variables at each node 20 (0-19)
        [varname, vartype, dimids, natts] = netcdf.inqVar(ncid,0);

        % load variables
        Temp = netcdf.getVar(ncid,6);       % temp in K
        Hum = netcdf.getVar(ncid,7);        % relative humidity in percent
        H = netcdf.getVar(ncid,8);          % Geopotential Height in m
        Plevs = netcdf.getVar(ncid,3);      % 37 pressure levels in Pa 
        netcdf.close(ncid)

        % flip the N-S axis
        Temp = flipdim(Temp,2);
        Temp = permute(Temp,[2,1,3]);
        Hum = flipdim(Hum,2);               % flip the N-S axis
        Hum = permute(Hum,[2,1,3]);    
        H = flipdim(H,2);                   % flip the N-S axis
        H = permute(H,[2,1,3]);

        % replacing the non-usefull values by NaN
        Temp(Temp==-999)=NaN;
        Hum(Hum==-999)=NaN;
        H(H==-999)=NaN;


        % getting the longitude and latitude coordinates
        numlat = size(Temp,1);
        numlon = size(Temp,2);
        ncid = netcdf.open(infile2,'NC_NOWRITE');
        lats = netcdf.getVar(ncid,1);
        latgrid = repmat(lats,[1,1,37]);
        latgrid = permute(latgrid,[2,1,3]);
        latgrid = flipdim(latgrid,1);
        lons = netcdf.getVar(ncid,2);
        longrid = repmat(lons,[1,1,37]);
        longrid = permute(longrid,[2,1,3]);
        longrid = flipdim(longrid,1);
        % close netcdf file
        netcdf.close(ncid)


        % Load pressure levels and replicate vector to make 3d matrix same size as variables
        Plevs = Plevs./100; % convert into hPa (or millibars) 1000 hPa at surface
        Pressure = repmat(Plevs,[1,numlat,numlon]);
        Pressure = permute(Pressure,[2,3,1]);


        % Get list of points to look at in analysis
        [xx,yy] = meshgrid(1:numlon,1:numlat);
        latlist = reshape(latgrid(:,:,1),[],1);
        lonlist = reshape(longrid(:,:,1),[],1);
        xlist = reshape(xx,[],1);
        ylist = reshape(yy,[],1);


 
        % When the landmask is coarser it can introduce extreme values at
        % sharp transitions. To cope with this I perform an interpoaltion
        % over the land points, such the land mask is removed.
        % interpolate the land 
        n_grid_points = size(Pressure,1)*size(Pressure,2);
        ix_vert =find(sum(sum(isnan(Temp)))==n_grid_points);
        if isempty(ix_vert)
            ix_vert = size(Temp,3);
        end
        longrid_vector = reshape(longrid(:,:,1),[],1);
        latgrid_vector= reshape(latgrid(:,:,1),[],1);
        
        % check for land mask only those layers that have land within them.
        parfor kkk=1:ix_vert(1)-1 
             ix = find(isnan(Temp(:,:,kkk))==1);
             if ~isempty(ix)
                 Temp_temp = reshape(Temp(:,:,kkk),[],1);
                 Hum_temp = reshape(Hum(:,:,kkk),[],1);
                 H_temp = reshape(H(:,:,kkk),[],1);

                 longrid_vector_temp = longrid_vector;
                 latgrid_vector_temp = latgrid_vector;
                 
                 % the land points that need to be interpolated
                 longrid_vector_needed = longrid_vector(ix);
                 latgrid_vector_needed = latgrid_vector(ix);
                 
                 % removing the land poitns with currently NaN values
                 Temp_temp(ix)=[];
                 Hum_temp(ix)=[];
                 H_temp(ix)=[];
                 
                 longrid_vector_temp(ix)=[];
                 latgrid_vector_temp(ix)=[];
                 

                 % do the interpolation for the landpoints
                 Temp_vector_needed = griddata(double(longrid_vector_temp),double(latgrid_vector_temp),double(Temp_temp),double(longrid_vector_needed),double(latgrid_vector_needed),'linear'); 
                 Hum_vector_needed = griddata(double(longrid_vector_temp),double(latgrid_vector_temp),double(Hum_temp),double(longrid_vector_needed),double(latgrid_vector_needed),'linear'); 
                 H_vector_needed = griddata(double(longrid_vector_temp),double(latgrid_vector_temp),double(H_temp),double(longrid_vector_needed),double(latgrid_vector_needed),'linear'); 
                 Temp_new = Temp(:,:,kkk);
                 Hum_new = Hum(:,:,kkk);  
                 H_new = H(:,:,kkk);       
                 
                 % replacing the values
                 Temp_new(ix)=Temp_vector_needed;
                 Hum_new(ix)=Hum_vector_needed;
                 H_new(ix)=H_vector_needed;
                 
                 % storing them back in the array
                 Temp(:,:,kkk) = Temp_new;
                 Hum(:,:,kkk) = Hum_new;
                 H(:,:,kkk) =  H_new;
             end
        end
        clear Temp_new Hum_new H_new ix Temp_temp Hum_temp H_temp longrid_vector_temp latgrid_vector_temp       
        
        % select the grid points within the area of interest
        % the projection of the grid is not uniform for the longitude as the center is defined for a specific longitude coordinate
        latlist_temp = sort(abs(diff(latlist)));
        lat_res = abs(latlist_temp(ceil(length(latlist_temp)*0.95)))*1.5;

        lonlist_temp = diff(lonlist);
        lon_res = max(abs(lonlist_temp))*1.5;


        % make sure the grid selected is covered and can be interpolated by WRF. Add a resolution cell to make sure it is the case.
        ix = find(ymin-lat_res<= latlist & latlist <= ymax+lat_res  & xmin-lon_res<= lonlist & lonlist<= xmax+lon_res) ;
        xlist = xlist(ix);
        ylist = ylist(ix);
        latlist = latlist(ix);
        lonlist = lonlist(ix);
        latlist=double(latlist);
        lonlist=double(lonlist);
        
        numy = length(unique(latlist));
        numx = length(unique(lonlist));
        ulatlist = unique(latlist);
        ulonlist = unique(lonlist);
        uxlist = unique(xlist);
        uylist = unique(ylist);
        
        

       
        
        if save_complete==0
            % saving the information for support plotting 
            wrf.wrf_lonlat = [lonlist latlist];
            wrf.region = [[xmin xmin xmax xmax xmin]' [ymin ymax ymax ymin ymin]'];
            deminfo.xmin = xmin;
            deminfo.xmax = xmax;
            deminfo.ymax = ymax;
            deminfo.ymin = ymin;
            deminfo.dem = dem;
            wrf.deminfo = deminfo;
            clear deminfo
            % checking if the file already exist. Yes append otherwiuze create it
            if exist('tca_support.mat','file')==2
              save('tca_support.mat','-append','wrf')
            else
              save('tca_support.mat','wrf')        
            end
            save_complete=1;
        end


               
        % Could not find the wrf used equation as they appear to be mixed
        % with latent heat etc. Istead I used the equations used in ERA-I
        % (see IFS documentation part 2: Data assimilation (CY25R1)). Calculate 
        % saturated water vapour pressure (svp) for water (svpw) using Buck 1881  
        % and for ice (swpi) from Alduchow and Eskridge (1996) euation AERKi
        svpw = 6.1121.*exp((17.502.*(Temp-273.16))./(240.97+Temp-273.16));
        svpi = 6.1121.*exp((22.587.*(Temp-273.16))./(273.86+Temp-273.16));
        tempbound1 = 273.16; %0
        tempbound2 = 250.16; %-23    
        

        svp = svpw;
        % Faster expression
        wgt = (Temp - tempbound2)/(tempbound1 - tempbound2);
        svp = svpi+ (svpw - svpi).*wgt.^2; % not sure why have wgt^2, Romain does it in his pyaps code.
        ix_bound1 = find(Temp > tempbound1);
        svp(ix_bound1) = svpw(ix_bound1);
        ix_bound2 = find(Temp < tempbound2);
        svp(ix_bound2) = svpi(ix_bound2);

        e = Hum./100.*svp;
        % NB, P = Pressure in hydrostatic equation, NOT dry partial pressure as
        % stated (incorrectly) in Jolivet et al. 2011.
        P = Pressure;


        %% Convert Geopotential to Geopotential Height and then to Geometric Height
        g0 = 9.80665;

        % map of g with latitude
        g = 9.80616.*(1 - 0.002637.*cosd(2.*latgrid) + 0.0000059.*(cosd(2.*latgrid)).^2);
        % map of Re with latitude
        Rmax = 6378137; 
        Rmin = 6356752;
        Re = sqrt(1./(((cosd(latgrid).^2)./Rmax^2) + ((sind(latgrid).^2)./Rmin^2)));        


        %% Calculate Geometric Height, Z
        Z = (H.*Re)./(g/g0.*Re - H);
        
        % Find middle of scene to work out glocal and Rlocal for use later
        midx = round(mean(uxlist));
        midy = round(mean(uylist));
        glocal = g(midy,midx,1);
        Rlocal = Re(midy,midx,1);
        
        

        %% Find middle of scene to work out glocal and Rlocal for use later
        cdslices = maxdem/vertres +1;
        cdstack = zeros(size(lonlist,1),cdslices);
        cdstack_dry = zeros(size(lonlist,1),cdslices);
        cdstack_wet = zeros(size(lonlist,1),cdslices);

        XI=(0:zincr:zref)';
        gh = glocal.*(Rlocal./(Rlocal+XI)).^2; %gravity with height for XI height range

        % Interpolate Temp P and e from 0:20:15000 m
        % then integrate using trapz to estimate delay as function of height
        parfor m = 1:size(lonlist,1);
            xn = xlist(m);
            yn = ylist(m);

            %interpolate at res zincr before doing integration
            X=double(squeeze(Z(yn,xn,:)));
            YP=double(squeeze(P(yn,xn,:)));
            Ye=double(squeeze(e(yn,xn,:)));
            YT=double(squeeze(Temp(yn,xn,:)));


            YT(isnan(Ye))=[];
            YP(isnan(Ye))=[];
            X(isnan(Ye))=[];
            Ye(isnan(Ye))=[];

            YeI = interp1(X,Ye,XI,'spline')*100; %change from hPa to Pa
            YPI = interp1(X,YP,XI,'spline')*100; %change from hPa to Pa
            YTI = interp1(X,YT,XI,'spline');
            

            tmp1 = ((k2-(Rd*k1/Rv)).*YeI./YTI + k3.*YeI./(YTI.^2));
            Lw = (10^-6).*-1*flipud(cumtrapz(flipud(XI),flipud(tmp1)));
            % L is zenith delay one way in metres
            tmp2 = k1.*YPI./YTI;
            %Ld = (10^-6).*-1*flipud(cumtrapz(flipud(XI),flipud(tmp2)));             % This is using P/T expression (Hanssen, 2001)
            gm = glocal; 
            Ld = (10^-6).*((k1*Rd/gm).*(YPI - YPI(zref/zincr +1)));                 % This is P0 expression (Hanssen, 2001)

            % Interpolate important part (i.e. total delay at elevations
            % less than maxdem) at high res i.e. vertres, and put in cdstack.
            cdI=(0:vertres:maxdem)';
            LdI=interp1(XI,Ld,cdI,'spline');
            LwI=interp1(XI,Lw,cdI,'spline');

            %cdstack(j,i,:) = LsI;
            cdstack_dry(m,:) = LdI;
            cdstack_wet(m,:) = LwI;
        end

        % Interpolate each cdstack layer onto a grid given by the DEM extents
        % in UTM m.
        xsmpres = (xmax-xmin)/nncols;
        ysmpres = (ymax-ymin)/nnrows;
        [xi,yi] = meshgrid(xmin+0.5*xsmpres:xsmpres:xmax-0.5*xsmpres,ymin+0.5*ysmpres:ysmpres:ymax-0.5*ysmpres);

        cdstack_interp_dry = zeros(nnrows,nncols,cdslices);
        parfor n = 1:cdslices
           slicelist = reshape(cdstack_dry(:,n),[],1);
           newslice = griddata(lonlist,latlist,slicelist,xi,yi,'linear'); 
           %F = TriScatteredInterp(lonlist,latlist,slicelist,'linear');
           %newslicelist = F(xi,yi);
           %newslice = reshape(newslicelist,nnrows,nncols);
           cdstack_interp_dry(:,:,n)= flipud(newslice); % flipud due to ypositive downpage for matlab figs, and ypositive uppage for utmy
        end    

        cdstack_interp_wet = zeros(nnrows,nncols,cdslices);
        parfor n = 1:cdslices
           slicelist = reshape(cdstack_wet(:,n),[],1);
           newslice = griddata(lonlist,latlist,slicelist,xi,yi,'linear'); 
           %F = TriScatteredInterp(lonlist,latlist,slicelist,'linear');
           %newslicelist = F(xi,yi);
           %newslice = reshape(newslicelist,nnrows,nncols);
           cdstack_interp_wet(:,:,n)= flipud(newslice); % flipud due to ypositive downpage for matlab figs, and ypositive uppage for utmy
        end  
        % keeping the coordinates in the same grid as the data
        xi = flipud(xi);
        yi = flipud(yi);


        % Pull out delays from cdstack layers that match dem heights
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
        
        
        % note that this is a one way Zenith delay and not a slant delay.
        % Units are in cm
        % Output wet correction
        wetcorrection = wetcorrection*100;                 % delay in [cm]
        fid = fopen(outfile,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(wetcorrection,[],1)]';
        tally = fwrite(fid,data_write,'double');
        fclose(fid);
        fid = fopen(outfile2,'w');
        fwrite(fid,wetcorrection','real*4');
        fclose(fid);   

        % Output dry correction
        hydrcorrection = hydrcorrection*100;                 % delay in [cm]  
        fid = fopen(hydroutfile,'w');
        data_write =  [reshape(xi,[],1) reshape(yi,[],1) reshape(hydrcorrection,[],1)]';
        tally = fwrite(fid,data_write,'double');
        fclose(fid);     
        fid = fopen(hydroutfile2,'w');
        fwrite(fid,hydrcorrection','real*4');
        fclose(fid);     

        
        fprintf([num2str(d) ' completed out of ' num2str(n_dates) '\n' ])
        
    else
        fprintf(['No wrf data for ' date ' \n'])
    end
end
