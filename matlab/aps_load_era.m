function [ Temp,WVapour,Geopot,Pressure,longrid,latgrid,xx,yy,lon0360_flag] = aps_load_era(file,era_data_type) 
% loading ERA-I or ERA5 data from ECMWF website or ERA-I from BADC website
% Bekaert David
% modifications
% DB    10/04/2016  extract code from aps_era_SAR.m to make code modular
% DB 	07/06/2017  Update syntax to include ERA5 model

%%% Example on how to load netcdf files 
% ncid = netcdf.open(file,'NC_NOWRITE');
% [numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);
% [dimname, dimlen] = netcdf.inqDim(ncid,0); 
%
%         for k=1:numvars
%             [dimname, dimlen] = netcdf.inqVar(ncid,k-1); 
%             fprintf([num2str(k-1) ' - ' dimname '\n'])
%         end

% debug figure to test and validate dataloading.
debug_fig = 0;

% open the netcdf
ncid = netcdf.open(file,'NC_NOWRITE');

% read netcdf variables and get number of variables
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);      


%% Swapping between BADC and ECMWF website data
if strcmpi(era_data_type,'ECMWF')
    % ECMWF data has field data and a scale plus offset.
    % Depending if this exist its added to the data
    for i = 0:numvars-1          
        [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
        flag = 0;
        for j = 0:numatts - 1
            attname1 = netcdf.inqAttName(ncid,i,j);
            attname2 = netcdf.getAtt(ncid,i,attname1);

            if strcmp('add_offset',attname1)
                offset = attname2;
            end

            if strcmp('scale_factor',attname1)
                scale = attname2;
                flag = 1;
            end
        end

        if flag
            eval([varname '= double(netcdf.getVar(ncid,i))*scale + offset;'])
        else
            eval([varname '= double(netcdf.getVar(ncid,i));'])
        end
        clear varname xtype dimids numatts scale offset      
    end


    % flip along third dimension, ECMWF format flipped compared to BADC 
    Temp = t(:,:,end:-1:1);       
    Hum = r(:,:,end:-1:1);
    Geopot = z(:,:,end:-1:1);
    Plevs = flipud(level);


elseif strcmpi(era_data_type,'BADC')
    % This is for ERA-I from BADC
    % Datafield are structured differently with different names
    % than ECMWF website data.

    % variables at each node 20 (0-19)
    [varname, vartype, dimids, natts] = netcdf.inqVar(ncid,0);

    % load data for variables 5 (temp), 11 (relative humidity) and 4 (geopotential height)
    Temp = double(netcdf.getVar(ncid,5));       % temp in K
    Hum = double(netcdf.getVar(ncid,11));       % relative humidity in percent
    Geopot = double(netcdf.getVar(ncid,4));     % Geopotential in m^2/s^2
    Plevs = double(netcdf.getVar(ncid,2));      % 37 pressure levels in hPa (or millibars) 1000 hPa at surface                   
end

% Permute to a 3 D grid
Temp = permute(Temp,[2,1,3]);       
Hum = permute(Hum,[2,1,3]);
Geopot = permute(Geopot,[2,1,3]);


% Same for lats and lons
lats = netcdf.getVar(ncid,1);
n_latitude_points = size(lats,1);
lons = netcdf.getVar(ncid,0);
n_longitude_points = size(lons,1);


% close nedtcdf
netcdf.close(ncid)


% adapting to the right lon lat grid size
Pressure = repmat(Plevs,[1,n_latitude_points,n_longitude_points]);
Pressure = permute(Pressure,[2,3,1]);


latgrid = repmat(lats,[1,37,n_longitude_points]);
latgrid = permute(latgrid,[1,3,2]);
longrid = repmat(lons,[1,37,n_latitude_points]);
longrid = permute(longrid,[3,1,2]);

% Get list of points to look at in analysis
[xx,yy] = meshgrid(1:n_longitude_points,1:n_latitude_points);


% (see IFS documentation part 2: Data assimilation (CY25R1)). 
% Calculate saturated water vapour pressure (svp) for water
% (svpw) using Buck 1881 and for ice (swpi) from Alduchow
% and Eskridge (1996) euation AERKi
svpw = 6.1121.*exp((17.502.*(Temp-273.16))./(240.97+Temp-273.16));
svpi = 6.1121.*exp((22.587.*(Temp-273.16))./(273.86+Temp-273.16));
tempbound1 = 273.16; %0
tempbound2 = 250.16; %-23      
svp = svpw;

% Faster expression
wgt = (Temp - tempbound2)/(tempbound1 - tempbound2);
svp = svpi+ (svpw - svpi).*wgt.^2;
ix_bound1 = find(Temp > tempbound1);
svp(ix_bound1) = svpw(ix_bound1);
ix_bound2 = find(Temp < tempbound2);
svp(ix_bound2) = svpi(ix_bound2);
WVapour = Hum./100.*svp;
clear Hum

% inform about the organisation of the longitudes
if sum(lons>180)>1
    lon0360_flag = 'y';
else
    lon0360_flag = 'n';
end

% validation plots
if debug_fig ==1
    figure('position',[ 3         628        1402         586]);
    subplot(2,5,1)
    imagesc(Temp(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('temp upper atmo')
    subplot(2,5,2)
    imagesc(Temp(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('temp lower atmo')

    subplot(2,5,3)
    imagesc(Pressure(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('pressure upper atmo')
    subplot(2,5,4)
    imagesc(Pressure(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('pressure lower atmo')



    subplot(2,5,6)
    imagesc(WVapour(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('Water vapour upper atmo')
    subplot(2,5,7)
    imagesc(WVapour(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('Water vapour lower atmo')


    subplot(2,5,8)
    imagesc(Geopot(:,:,end))
    colorbar
    axis xy
    axis equal
    axis tight
    title('geopotential upper atmo')
    subplot(2,5,9)
    imagesc(Geopot(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('geopotential lower atmo')


    subplot(2,5,5)
    imagesc(latgrid(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('lat')
    subplot(2,5,10)
    imagesc(longrid(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('lon')
end



