function [] = aps_spectrometer_PWV_meris(batchfile)
% aps_spectrometer_PWV_comparison(datelist)
% Scipt to load meris data, mask out clouds. The meris data is assumed to 
% be structured in date folders. The batchfile contains the full path to the
% meris files in these folders. Note that the first line of the batchfile should read "files".
%
%
% INPUTS:
% batchfile             A txt file containing the full path and file names of the
%                       meris data that needs to be processed. The first
%                       line of this file should read "files". The data
%                       should be structured in date folders.
% xlims                 Limits in the x-direction, either in degrees
% ylims                 Limits in the y-direction, either in degrees
% smpres                The output resolution, either in degrees
%                       Units needs to be consistend with xlims and ylims.
%
%
% By David Bekaert - University of Leeds
% August 2014

fig_test = 1;                   % when 1 show the dem as debug figure



if nargin<1
    fprintf('aps_spectrometer_PWV_meris(batchfile) \n')
    error('myApp:argChk', ['Not enough input arguments...\n'])
end


smpres = getparm_aps('region_res',1); % in degrees
xlims = getparm_aps('region_lon_range',1);
ylims = getparm_aps('region_lat_range',1);
stamps_processed = getparm_aps('stamps_processed',1);



%% the actual scripting
%bounds for all ifgs in degrees or in meters
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);


% getting the number of files to be processed
files = char(textread(batchfile,'%s','headerlines',1));
ndates = size(files,1);

% loading the date information  
if strcmp(stamps_processed,'y')
   ps = load(getparm_aps('ll_matfile',1));
   ifgs_dates = ps.day;
   fprintf('Stamps processed structure \n')
else
    ifgday_matfile = getparm_aps('ifgday_matfile',1);
    ifgs_dates = load(ifgday_matfile);
    ifgs_dates = ifgs_dates.ifgday;
    ifgs_dates = reshape(ifgs_dates,[],1);
    ifgs_dates = unique(ifgs_dates);
end


fprintf(['Lon range: ' num2str(xmin) ' -- '  num2str(xmax) ' degrees\n'])
fprintf(['Lat range: ' num2str(ymin) ' -- '  num2str(ymax) ' degrees\n'])
fprintf(['Output resolution is assumed to be ' num2str(smpres) ' degrees \n'])


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



%start loop here to calculate atmos correction for each date
fprintf('Starting the masking and writing of the PWV for each SAR date \n')
for n = 1:ndates
    
    file = [files(n,:)];
    outfile_watervapor = [pathlist{n} filesep datelist(n,:) '_ZPWV_nointerp.xyz'];
    
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


    %% Calculate wet delay if desired

    %do masking, interpolating etc here.
    %convert flag matrix to long list of binary numbers e.g. 1000010000010000
    flagbin = dec2bin(meris(:,:,33));


    %% Generating cloud mask
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

    %% loading and masking of the data
    %load water vapour
    watervap=meris(:,:,14);   % g/cm^2  water vapour
    %mask out dodgy pixels
    corwatervap = watervap .*maskmat;     

    %% mask water additional
    %%% aditonal masking
    % water mask
    water = bin2dec(flagbin(:,3));
    water=reshape(water,mer_rows,mer_cols);

    %decrease water mask edges, to make sure we do not remove land points
    se = strel('square', 3);
    water_new = imerode(water,se);
    water_new(water_new==1)=NaN;
    water_new(water_new==0)=1;
    clear water
    corwatervap = corwatervap.*water_new; 



    %% saving the data
    fid = fopen('out.bin','w');
    corwatervap = corwatervap';   
    fwrite(fid,corwatervap,'real*4');
    fclose(fid);

    % do gmt resample and cut
    xyz2grd_cmd = ['xyz2grd -R',num2str(mer_xmin),'/',num2str(mer_xmax),'/',num2str(mer_ymin),'/',num2str(mer_ymax),' -I',num2str(mer_cols),'+/',num2str(mer_rows),'+ out.bin -Gtmp.grd  -F -ZTLf'];
    [a,b] = system(xyz2grd_cmd);

    % downsample and output without interpolation files
    grdsmp_cmd = ['grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' -F tmp.grd -Gtmp_smp.grd'];
    [a,b]=system(grdsmp_cmd);
    grd2xyz_cmd = ['grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo >' outfile_watervapor];
    [a,b]=system(grd2xyz_cmd);




    % load gaussian and nointerp back in, combine, and output as outfile_gauss
    % opening the data file (not-interpolated)
    nointfid = fopen(outfile_watervapor,'r');
    data_vector = fread(nointfid,'double');
    fclose(nointfid);
    % reshaping into the right n column matrix
    data = reshape(data_vector,3,[])';
    noint = data(:,3);
    xy = data(:,[1:2]);
    clear data data_vector

%     figure; scatter3(xy(:,1),xy(:,2),noint,15,noint,'filled'); view(0,90); axis equal; axis tight
    
    % writing out the date again as a binary table
    data_write = [xy noint]';
    clear noint xy
    %output
    fid = fopen(outfile_watervapor,'w');
    fwrite(fid,data_write,'double');
    fclose(fid);
    clear data_write
    
    fprintf([num2str(n) ' completed out of ' num2str(ndates) '\n'])
    
end

[a b] = system('!rm tmp.grd tmp.xyz out.bin tmp2.grd tmp_fil.grd tmp_smp.grd tmp_smp.xyz tmp_smp2.grd tmp_smp_fil.grd');

