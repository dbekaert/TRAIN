function [] = aps_meris_SAR(batchfile)
% aps_meris_SAR(datelist)
% Scipt to load meris data, mask out clouds, interpolate over gaps and cut 
% the data to the right size for which the tropospheric delay is being computed. 
% The DEM file inputed should have an asociated
% ".rsc" file, with the same filename as the DEM. The ".rsc" files should
% contain a WIDTH, LENGTH, X_FIRST, Y_FIRST, X_STEP, Y_STEP and optional a 
% FORMAT string. The meris data is assumed to be structured in date folders. 
% The batchfile contains the full path to the meris files in these folders. 
% Note that the first line of the batchfile should read "files".
%
%
% INPUTS:
% batchfile             A txt file containing the full path and file names of the
%                       meris data that needs to be processed. The first
%                       line of this file should read "files". The data
%                       should be structured in date folders.
% demfile               Full path to the DEM file. The DEM needs to me in meters.
% xlims                 Limits in the x-direction, either in degrees
% ylims                 Limits in the y-direction, either in degrees
% wetdry                1=calc wet only, 2=calc dry only, 3=calc both (default).
% demnull               The value for no DEM data, default is -32768.
% smpres                The output resolution, either in degrees
%                       Units needs to be consistend with xlims and ylims.
%
% OPTIONAL INPUTS:
% conversion            Conversion factor, default 6.2 .
% scaleheight           Scale height, default 8340 m.
%
% OUTPUTS:
% Depending on your selected option it will give computed dry, wet, or
% combined delay maps in cm. By default the wet delay is
% computed as the dry delay is not recomended to be used for meris, use era
% or weather model instead
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
% Modified from Richard Walters - Oxford / Leeds University 2012-2013
% Modifications:
% DB 	02/2013		Convert script to a function and add syntax
% DB    02/2013		Allow for a more flexible input file name
% DB    02/2013     Make the output files independant of grd files.
% DB    02/2013     Allow the precision of the DEM to be specified in the
%                   dem input file (.rsc file).
% DB    02/2013     Adding extra syntax to the code
% DB    04/2013     Adding extra check for input arguments
% DB    04/2013     Incorporate getparm_aps from the aps_toolbox.
% DB    10/2013     Changed filename to be more consistent with toolbox
% DB    03/2014     Suppress command line output
% DB    03/2014     Include an extra check for the DEM grd file and the
%                   selected crop
% DB    07/2014     Include user interaction message for DEM check
% DB    07/2014     Redefine meris_lat(lon)_range to region_lat(lon)_range
% DB    08/2014     Include option to have factors varying for each SAR date
% DB    02/2015     Remove the hydrostatic component and only keep wet delay.
% DB    03/2015     Remove scale height paramter

                                
if nargin<1
    fprintf('load_meris(batchfile) \n')
    error('myApp:argChk', ['Not enough input arguments...\n'])
end


conversion_vector = getparm_aps('spectrometer_PIconversion');
smpres = getparm_aps('region_res'); % in degrees

xlims = getparm_aps('region_lon_range');
ylims = getparm_aps('region_lat_range');
dryconversion = 0.23;
wetdry = 1;                         % 1=calc wet only, 2=calc dry only, 3=calc both (default)



%% the actual scripting
%bounds for all ifgms in degrees or in meters
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);
fprintf(['Lon range: ' num2str(xmin) ' -- '  num2str(xmax) ' degrees\n'])
fprintf(['Lat range: ' num2str(ymin) ' -- '  num2str(ymax) ' degrees\n'])
fprintf(['Output resolution is assumed to be ' num2str(smpres) ' degrees \n'])


% getting the number of files to be processed
files = char(textread(batchfile,'%s','headerlines',1));
ndates = size(files,1);

% loading the date information  
stamps_processed = getparm_aps('stamps_processed');
if strcmp(stamps_processed,'y')
   ps = load(getparm_aps('ll_matfile'));
   ifgs_dates = ps.day;
   fprintf('Stamps processed structure \n')
else
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    ifgs_dates = reshape(ifgs_dates,[],1);
    ifgs_dates = unique(ifgs_dates);
end


% extracting the dates from the filenames
for k=1:ndates
    [path,filename_temp,ext_temp] = fileparts(files(k,:));
    clear filename_temp ext_temp
    [path_temp,date,ext_temp] = fileparts(path);
    clear path_temp ext_temp
    
    % save the paths as structures to allow for variable path lengths
    pathlist{k} = path;
    clear path
    % saving the date information
    datelist(k,:) =date;
    clear date
end


if length(conversion_vector)>1
   fprintf('Conversion factor varies for each SAR date.\n')
end


%start loop here to calculate atmos correction for each date
fprintf('Starting the computation for each SAR date \n')
for n = 1:ndates
    if length(conversion_vector)>1
        ix_date_postion = find(ifgs_dates==datenum(datelist(n,:),'yyyymmdd'));
        conversion =conversion_vector(ix_date_postion);
    else
       conversion =conversion_vector;
    end
    
    file = [files(n,:)];
    outfile = [pathlist{n} filesep datelist(n,:) '_SWD_nointerp.xyz'];
    outfile_gauss = [pathlist{n} filesep datelist(n,:) '_SWD_gauss.xyz'];
    outfile_surf = [pathlist{n} filesep datelist(n,:) '_SWD_surf.xyz'];
    dryoutfile = [pathlist{n} filesep datelist(n,:) '_SHD.xyz'];
    
    % get info on meris file from tif
    tifinfo = geotiffinfo(file);
    mer_xmin = tifinfo.CornerCoords.X(1);
    mer_xmax = tifinfo.CornerCoords.X(2);
    mer_ymin = tifinfo.CornerCoords.Y(3);
    mer_ymax = tifinfo.CornerCoords.Y(1);
    mer_cols = tifinfo.Width;
    mer_rows = tifinfo.Height;
    mer_x = (mer_xmax - mer_xmin)/mer_cols;
    mer_y = (mer_ymax - mer_ymin)/mer_rows;
    % read in tif
    meris=imread(file,'tif');


    %% Calculate wet delay
    %do masking, interpolating etc here.
    %convert flag matrix to long list of binary numbers e.g. 1000010000010000
    flagbin = dec2bin(meris(:,:,33));

    %take certain bits of binary numbers (e.g. second bit, 1=cloud), and reshape
    %into matrix mask
    cloud = bin2dec(flagbin(:,2));
    cloudmat=reshape(cloud,mer_rows,mer_cols);
    cloudmat(cloudmat==1)=NaN;

    % turn ones (i.e. mask pixels) into NaNs so when we multiply masks together,
    % NaNs will penetrate through the stack
    pconf = bin2dec(flagbin(:,23));
    pconfmat=reshape(pconf,mer_rows,mer_cols);
    pconfmat(pconfmat==1)=NaN;

    lowp = bin2dec(flagbin(:,24));
    lowpmat=reshape(lowp,mer_rows,mer_cols);
    lowpmat(lowpmat==1)=NaN;


    maskmat= cloudmat .*pconfmat .*lowpmat;
    maskmat(isnan(maskmat))=1;

    %extend mask edges
    se = strel('square', 3);
    maskmat = imdilate(maskmat,se);

    maskmat(maskmat==1)=NaN;
    maskmat(maskmat==0)=1;

    %load water vapour
    watervap=meris(:,:,14);
    %convert from zenith to slant using view_zenith tie point grid
    watervap=watervap./cosd(meris(:,:,42));
    %convert from g/cm^2 slant water vapour to cm phase delay 
    watervap=watervap.*conversion;
    %mask out dodgy pixels
    corwatervap = watervap .*maskmat;     


%         %%% aditonal masking
%         % water mask
%         water = bin2dec(flagbin(:,3));
%         water=reshape(water,mer_rows,mer_cols);
%         
%         %decrease water mask edges, to make sure we do not remove land points
%         se = strel('square', 3);
%         water_new = imerode(water,se);
%         water_new(water_new==1)=NaN;
%         water_new(water_new==0)=1;
%         clear water
%         
%         corwatervap = corwatervap.*water_new; 
%         
%         
%         %%%% mask for outliers in the estimates
%         temp_data = (corwatervap - medfilt2(corwatervap, [3 3]));
%         temp  = reshape(temp_data,[],1);
%         temp(isnan(temp))=[];
%         outlier_mask = abs(temp_data)>4*std(temp);
%         outlier_mask = outlier_mask+1;
%         outlier_mask(outlier_mask==2)=NaN;
%         corwatervap = corwatervap.*outlier_mask; 
%         clear outlier_mask water maskmat


    fid = fopen('out.bin','w');
    corwatervap = corwatervap';   
    fwrite(fid,corwatervap,'real*4');
    fclose(fid);


    % do gmt resample and cut
    xyz2grd_cmd = ['xyz2grd -R',num2str(mer_xmin),'/',num2str(mer_xmax),'/',num2str(mer_ymin),'/',num2str(mer_ymax),' -I',num2str(mer_cols),'+/',num2str(mer_rows),'+ out.bin -Gtmp.grd  -F -ZTLf'];
    [a,b] = system(xyz2grd_cmd);


    % export tmp.grd as tmp.xyz
    grd2xyz_cmd = ['grd2xyz -R',num2str(mer_xmin),'/',num2str(mer_xmax),'/',num2str(mer_ymin),'/',num2str(mer_ymax),' tmp.grd -bo > tmp.xyz'];
    [a,b]=system(grd2xyz_cmd);


    % gaussian filter tmp.grd to get tmp2.grd
    grdfil1_cmd = 'grdfilter tmp.grd -Gtmp2.grd -Ni -D2 -Fg50'; %10000
    [a,b]=system(grdfil1_cmd);


    % use surface to interpolate using tmp.xyz and create tmp_fil.grd
    grdfil_cmd = ['surface -R',num2str(mer_xmin),'/',num2str(mer_xmax),'/',num2str(mer_ymin),'/',num2str(mer_ymax),' tmp.xyz -I',num2str(mer_cols),'+/',num2str(mer_rows), '+ -bi -Gtmp_fil.grd -T0.5'];
    [a,b]=system(grdfil_cmd);

    % downsample and output nointerp, gaussian and surface files
    %nointerp
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp.grd -Gtmp_smp.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo >' outfile];
    [a,b]=system(grd2xyz_cmd);

    %gauss
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp2.grd -Gtmp_smp2.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp2.grd -bo >' outfile_gauss];
    [a,b]=system(grd2xyz_cmd);
    %surf
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp_fil.grd -Gtmp_smp_fil.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp_fil.grd -bo >' outfile_surf];
    [a,b]=system(grd2xyz_cmd);

    % load gaussian and nointerp back in, combine, and output as outfile_gauss
    % opening the data file (not-interpolated)
    nointfid = fopen(outfile,'r');
    data_vector = fread(nointfid,'double');
    fclose(nointfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    noint = data(:,3);
    clear data data_vector


    % opening the gaussian interpolated file
    gaussfid = fopen(outfile_gauss,'r');
    data_vector = fread(gaussfid,'double');
    fclose(gaussfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    gauss = data(:,3);
    xy = data(:,[1:2]);
    clear data data_vector

    % (gaussian interpolated points only where no data available)
    noint(isnan(noint))=0;
    gauss(noint~=0)=0;
    gauss_interp=gauss+noint;
    gauss_interp(gauss_interp==0)=NaN;
    % writing out the date again as a binary table
    data_write = [xy gauss_interp]';
    clear gauss_interp noint gauss
    %output
    fid = fopen(outfile_gauss,'w');
    fwrite(fid,data_write,'double');
    fclose(fid);
    clear data_write




    


    %output incidence.xyz one time only
    if n==1
        incidence = meris(:,:,42)';
        fid = fopen('out.bin','w');
        fwrite(fid,incidence,'real*4');
        fclose(fid);

        xyz2grd_cmd = ['xyz2grd -R',num2str(mer_xmin),'/',num2str(mer_xmax),'/',num2str(mer_ymin),'/',num2str(mer_ymax),' -I',num2str(mer_cols),'+/',num2str(mer_rows),'+ out.bin -Gtmp.grd  -F -ZTLf'];
        [a,b]=system(xyz2grd_cmd);
        grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp.grd -Gtmp_smp.grd'];
        [a,b]=system(grdsmp_cmd);
        grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo > incidence.xyz'];
        [a,b]=system(grd2xyz_cmd);
    end
    
    fprintf([num2str(n) ' completed out of ' num2str(ndates) '\n'])
    
end

[a b] = system('!rm tmp.grd tmp.xyz out.bin tmp2.grd tmp_fil.grd tmp_smp.grd tmp_smp.xyz tmp_smp2.grd tmp_smp_fil.grd');

