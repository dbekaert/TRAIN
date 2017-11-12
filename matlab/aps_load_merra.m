function [ Temp,WVapour,Geopot,Pressure,longrid,latgrid,xx,yy,lon0360_flag] = aps_load_merra(file,merra_model,coeff) 
% loading merra data and output the variables for TRAIN
%
%     Copyright (C) 2016  Bekaert David
%     davidbekaert.com
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
% modifications
% DB    10/04/2016      Base this code on aps_load_era.m modular approach  
% DB    01/05/2016      Include MERRA2 support
% DB    12/11/2017      Change MERRA to be nc4 


% debug figure to test and validate dataloading.
debug_fig = 0;

% missing values in the hdf files of MERRA
missing_value = 999999986991104;

% specify the coeficients in case they are not given
if nargin<3
    coeff.Rv= 461.495;           % [j/kg/K] water vapour
    coeff.Rd= 287.05;            % [j/kg/K] specific gas constants dry air
end
if ~isfield(coeff,'Rd') || ~isfield(coeff,'Rv')
   error('Coefficients Rd and Rv are not spcified') 
end

%%
if strcmpi(merra_model,'merra')
    
    [a,b,file_ext] = fileparts(file);
    
    if strcmpi(file_ext,'.nc4')
        % open the netcdf
        ncid = netcdf.open(file,'NC_NOWRITE');

        % read netcdf variables and get number of variables
        [numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);      


         % loading all the variables from the netcdf
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
                eval([varname '= double(netcdf.getVar(ncid,i))*scale + offset;']);
            else
                eval([varname '= double(netcdf.getVar(ncid,i));']);
            end
            clear varname xtype dimids numatts scale offset      
         end

         % close nedtcdf
        netcdf.close(ncid)

        % definition of QV
        qv = QV;

        % definition of geopotential
        % Convert Geopotential Height to geopotential
        % important this is needed to be consitent with the other weather modules
        g0 = 9.80665;
        Geopot = H.*g0;

        % define geocoordinates
        lons = XDim;
        lats = YDim;
        clear XDim YDim 
        % define temperature
        Temp = T;
        clear T
        % define pressure
        Plevs = Height;
        clear Height

        % this is the time-stamps
        TIME./60;

    elseif  strcmpi(file_ext,'.hdf')
        %% OLD HDF5 files, seems no longer to work for new MERRA website
        % keep it for being backward compatible with earlier TRAIN version
        hdf_file = hdfinfo(file);
        Temp = hdfread(file,'/t');                      % temperature in [K]
        qv  = hdfread(file,'/qv');                      % specific humidity
        H =  hdfread(file,'/h');
        % Convert Geopotential Height to geopotential
        % important this is needed to be consitent with the other weather modules
        g0 = 9.80665;
        Geopot = H.*g0;                                 % geopotential
        Psurface = hdfread(file,'/ps');                 % surface pressure in [hPa]
        Plevs = hdfread(file,'/levels')';               % Pressure in [hPa]
        lons = hdfread(file,'/longitude')';
        lats = hdfread(file,'/latitude')';
    end
elseif strcmpi(merra_model,'merra2')
    
    % open the netcdf
    ncid = netcdf.open(file,'NC_NOWRITE');

    % read netcdf variables and get number of variables
    [numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);      


     % loading all the variables from the netcdf
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
            eval([varname '= double(netcdf.getVar(ncid,i))*scale + offset;']);
        else
            eval([varname '= double(netcdf.getVar(ncid,i));']);
        end
        clear varname xtype dimids numatts scale offset      
     end
     
     % close nedtcdf
    netcdf.close(ncid)
     
    % definition of QV
    qv = QV;
    
    % definition of geopotential
    % Convert Geopotential Height to geopotential
    % important this is needed to be consitent with the other weather modules
    g0 = 9.80665;
    Geopot = H.*g0;
    
    % define geocoordinates
    lons = lon;
    lats = lat;
    clear lon lat 
    % define temperature
    Temp = T;
    clear T
    % define pressure
    Plevs = lev;
    clear lev
       
    % this is the time-stamps
    time./60;
end
%%
% update the no data values with NaN's
Temp(Temp==missing_value)=NaN;
qv(qv==missing_value)=NaN;
Geopot(Geopot==missing_value)=NaN;


% number of lon lat grid nodes
n_latitude_points = size(lats,1);
n_longitude_points = size(lons,1);

% generate the pressure levels
Pressure = repmat(Plevs,[1,n_latitude_points,n_longitude_points]);
Pressure = permute(Pressure,[2,3,1]);


% change the order of variables such its consistent with the pressure levels
if strcmpi(merra_model,'merra2')
    Temp = permute(Temp,[2,1,3]);
    qv = permute(qv,[2,1,3]);
    Geopot = permute(Geopot,[2,1,3]);

elseif  strcmpi(merra_model,'merra')
    if strcmpi(file_ext,'.nc4')
        Temp = permute(Temp,[2,1,3]);
        qv = permute(qv,[2,1,3]);
        Geopot = permute(Geopot,[2,1,3]);
    elseif strcmpi(file_ext,'.hdf')
        Temp = permute(Temp,[2,3,1]);
        qv =  permute(qv,[2,3,1]);
        Geopot = permute(Geopot,[2,3,1]);
    end
end


% computation of the saturated water vapour e
% use mixing ratio r = E*e/(p-e), where E = Rd/Rv
% and specific humidity qv = r/(1+r)
E = coeff.Rd./coeff.Rv;
WVapour= qv.*Pressure./(E.*(1-qv)+qv);
clear qv

% Get list of points to look at in analysis
[xx,yy] = meshgrid(1:n_longitude_points,1:n_latitude_points);

% generate the lon lat grid
latgrid = repmat(lats,[1,length(Plevs),n_longitude_points]);
latgrid = permute(latgrid,[1,3,2]);
longrid = repmat(lons,[1,length(Plevs),n_latitude_points]);
longrid = permute(longrid,[3,1,2]);


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
    title('water vapour upper atmo')
    subplot(2,5,7)
    imagesc(WVapour(:,:,1))
    colorbar
    axis xy
    axis equal
    axis tight
    title('water vapour lower atmo')


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


% 
% 
% 
%  ncid = netcdf.open(file,'NC_NOWRITE');
% [numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);
% [dimname, dimlen] = netcdf.inqDim(ncid,0); 
% 
%         for k=1:numvars
%             [dimname, dimlen] = netcdf.inqVar(ncid,k-1); 
%             fprintf([num2str(k-1) ' - ' dimname '\n'])
%         end
