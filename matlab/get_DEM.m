function [dem,xmin,xmax,ymin,ymax,smpres,nncols,nnrows] = get_DEM()
% function to load the DEM, convert to a grid, resample to the user defined
% grid. Gives as output the DEM grid, the coordinates information, number
% of rows and columns of the grid.
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
% By Bekaert David - Unversity of Leeds
% modifications:
% 10/2015   DB  Adding support for grd file DEM's
% 10/2015   DB  bug fix for grd file, when the longitude is defined +360degree
% 10/2015   DB  Checking for GMT5 compatibility. No issues found, but
%               information is passed on
% 11/2015   DB  Add system call check for errors
% 02/2016   DB  Add a fix in case format is not given
%               Attemt to autofix between int16 and short format
% 04/2016   DB  Check if the DEM exist
% 04/2016   DB  Add option to give DEM with .xml file too
% 11/2017   DB  finnish xml implementation, allow for gmt string to be
%               given to call gmt outside matlab



% retreiving which version of GMT this is
[gmt5_above, gmt_version,GMT_string] = get_gmt_version;



fig_test = 1;               % when 1 show the dem as debug figure
smpdem = 'dem_smp.xyz';     % output file for the DEM

% information read from teh paramter file
demfile = getparm_aps('demfile');
pixelreg = 'node';%getparm_aps('pixelreg');      %added by hua wang, pixel registration mode (node or gridline)
% outfmt = getparm_aps('outfmt');
demnull = getparm_aps('dem_null');
xlims = getparm_aps('region_lon_range',1);
ylims = getparm_aps('region_lat_range',1);
smpres = getparm_aps('region_res',1);



%% no changes required below
% check if DEM exist
if exist(demfile,'file')~=2
    error([demfile ' does not exist'])
end

% check if the longitude and latitude are given correctly
xlims_test = [min(xlims) max(xlims)];
ylims_test = [min(ylims) max(ylims)];
if (sum([ylims_test-ylims])+sum([xlims_test-xlims]))>0
    error(['You region_lon_range and/or region_lat_range is not defined as [smallest largest], correct this.'])
end

% the region the delay needs to cover.
% slightly larger than the InSAR region
xmin = xlims(1);
xmax = xlims(2);
ymin = ylims(1);
ymax = ylims(2);

% getting the DEM path
[path_dem,filename_dem,ext_dem] = fileparts(demfile);

% checking if file exists
if exist(demfile,'file')~=2
    error('dem file does not exist')
end

% check if this is a grd file already
if strcmpi(ext_dem,'.grd')
    fprintf(['Specified DEM is a .grd file\n'])
    if exist('tmp.grd','file')==2
        delete('tmp.grd')
    end
    system_command = ['ln -s ' demfile ' tmp.grd'];
    aps_systemcall(system_command);
else

    clear filename_dem ext_dem
    file_bad=2;
    dem_support = [];
    if exist([demfile, '.rsc'],'file')==2
        file_bad=file_bad-1;
        dem_support = ['.rsc'];
    end
    if exist([demfile, '.xml'],'file')==2
        file_bad=file_bad-1;
        dem_support = ['.xml'];
    end
        
    if file_bad ==2
        fprintf(['Your DEM is not .grd.\n This is fine but you need a  .rsc or.xml file\n'])
        error([demfile   '.rsc or.xml not found' ])
    elseif file_bad ==0
        fprintf(['Found both ' demfile '.xml and .rsc. Will use .rsc file by default \n']);
        dem_support = ['.rsc'];
    end

    
    % gettign the DEM specifics
    if strcmpi(dem_support,'.xml')
        % get dem details from xml file
        [ncols] = get_parm_xml([demfile, '.xml'],'width');
        [nrows] = get_parm_xml([demfile, '.xml'],'length');
        [scheme] = get_parm_xml([demfile, '.xml'],'scheme');
        [number_bands] = get_parm_xml([demfile, '.xml'],'number_bands');
        [format] = get_parm_xml([demfile, '.xml'],'data_type');
        [delta_list,fields_delta] = get_parm_xml([demfile, '.xml'],'delta');
        [startingvalue_list,fields_startvalue] = get_parm_xml([demfile, '.xml'],'startingvalue');
        [endingvalue_list,fields_endvalue] = get_parm_xml([demfile, '.xml'],'endingvalue');
        
        
        xfirst = startingvalue_list(1);
        yfirst = startingvalue_list(2);
        xstep = delta_list(1);
        ystep = delta_list(2);
        
        
        % santiy check on the ending valuu
        xlast = xfirst + ncols*xstep;
        ylast = yfirst + nrows*ystep;
        if sum(abs(endingvalue_list-[xlast; ylast]))~0
            fprintf(['Did not expect the x-ylast value to be different from endingvalue_list, please check\n'])
            keyboard
        end
        
        if number_bands>1
            error('Only 1-band files are supported for DEM')
        end
        
        
    elseif strcmpi(dem_support,'.rsc')
        % get dem details from rsc file
        ncols_cmd = ['echo `grep WIDTH ', demfile, '.rsc | awk ''{print $2}''`>', path_dem ,filesep ,'temp'];
        aps_systemcall(ncols_cmd);
        nrows_cmd = ['echo `grep LENGTH ', demfile, '.rsc | awk ''{print $2}''`>>', path_dem ,filesep ,'temp'];
        aps_systemcall(nrows_cmd);
        xfirst_cmd = ['echo `grep X_FIRST ', demfile, '.rsc | awk ''{print $2}''`>>', path_dem ,filesep ,'temp'];
        aps_systemcall(xfirst_cmd);
        yfirst_cmd = ['echo `grep Y_FIRST ', demfile, '.rsc | awk ''{print $2}''`>>', path_dem ,filesep ,'temp'];
        aps_systemcall(yfirst_cmd);
        xstep_cmd = ['echo `grep X_STEP ', demfile, '.rsc | awk ''{print $2}''`>>', path_dem ,filesep ,'temp'];
        aps_systemcall(xstep_cmd);
        ystep_cmd = ['echo `grep Y_STEP ', demfile, '.rsc | awk ''{print $2}''`>>', path_dem ,filesep ,'temp'];
        aps_systemcall(ystep_cmd);
        format_cmd = ['echo `grep FORMAT ', demfile, '.rsc | awk ''{print $2}''`>', path_dem ,filesep ,'temp2'];
        aps_systemcall(format_cmd);

        % loading the DEM info data (numeric values)
        DEM_info = load([path_dem ,filesep 'temp']);
        aps_systemcall(['rm ' path_dem ,filesep 'temp']);
        ncols = DEM_info(1);
        nrows = DEM_info(2);
        xfirst = DEM_info(3);
        yfirst = DEM_info(4);
        xstep = DEM_info(5);
        ystep = DEM_info(6);
        clear DEM_info
        
        
        % loading the DEM info data (strings)
        fid = fopen([path_dem ,filesep 'temp2']) ;
        DEM_info2 = textscan(fid,'%s');
        fclose(fid);
        aps_systemcall(['rm ' path_dem ,filesep 'temp2']);
        format =  DEM_info2{1};
        clear DEM_info2
    end
    

    %% Getting the DEM precision
    %automatically checking dem format, by hua wang, 26 Feb 2015
    fid=fopen(demfile,'r');
    fseek(fid, 0, 'eof');
    pos=ftell(fid);
    byte=pos/ncols/nrows;
    fseek(fid, 0, 'bof');
    switch byte
      case 8 
          format_dem_str='d';     %double
      case 4 
          format_dem_str='f';     %float, real4
      case 2 
          format_dem_str='i';     %int16
      case 1 
          format_dem_str='h';     %int8
      otherwise
          error('no such dem format');
    end
    % give warning as we cannot discriminate between int16 and short
    if strcmpi(format_dem_str,'i')
        fprintf('Could be int16 or short format \n')
    end
    fclose(fid);

    

    format_dem_str_given = [];
    if ~isempty(format)
        if  strcmpi(format,'r4') || strcmpi(format,'real4')
           format_dem_str_given='f'; 
        elseif strcmpi(format,'h') || strcmpi(format,'short')
            format_dem_str_given='h';   
            % if the given fortmat is short then update the detection to be short
            if strcmpi(format_dem_str,'i')
                format_dem_str = 'h';
            end
        end
    end

    
    % checking if the specified precision is different, then that of the automated estimation
    if ~isempty(format_dem_str_given)
        if ~strcmpi(format_dem_str_given,format_dem_str)
            fprintf('The automated detected DEM precision is different then what you specified in dem.rsc file \n:')
            fprintf(['Yours: ' format_dem_str_given '\n'])
            fprintf(['Auto: ' format_dem_str '\n'])

            repeat=1;
            while repeat==1
                action_flag= str2num(input('Keep yours (1), Keep auto (2), Different (3)? [1, 2, or 3] ','s'));
                if isnumeric(action_flag)
                    if action_flag==1
                        format_dem_str = format_dem_str_given;
                        repeat=0;
                    elseif action_flag==2
                        format_dem_str = format_dem_str;
                        repeat=0;
                    elseif action_flag==3
                        action_flag= input('To what do you want to update this? [give a format recognised by GMT]','s');
                        format_dem_str = action_flag;
                        repeat=0;
                    end
                end
            end
            fprintf(['Please verify, and re-run'])
        end
    end


    %% getting the DEM as a grid again

    %getting the extends of the DEM
    %revised by hua wang, 26, Feb
    if strcmp(pixelreg,'gridline')==0
      xlast = xfirst + ncols*xstep;
      ylast = yfirst + nrows*ystep;
    else
      xlast = xfirst + (ncols-1)*xstep;
      ylast = yfirst + (nrows-1)*ystep;
    end

    if strcmp(pixelreg,'gridline')==0
      xyz2grd_cmd = [GMT_string 'xyz2grd -R',num2str(xfirst),'/',num2str(xlast),'/',num2str(ylast),'/',num2str(yfirst),' -I',num2str(ncols),'+/',num2str(nrows),'+ ',demfile,' -Gtmp.grd -N',num2str(demnull),' -F -ZTL' format_dem_str];
    else
      xyz2grd_cmd = [GMT_string 'xyz2grd -R',num2str(xfirst),'/',num2str(xlast),'/',num2str(ylast),'/',num2str(yfirst),' -I',num2str(ncols),'+/',num2str(nrows),'+ ',demfile,' -Gtmp.grd -N',num2str(demnull),' -ZTL' format_dem_str];
    end
    
    % down-sample dem to the grid as secified by the user with the given resolution    
    % in case the format is int16 lets first try to see if its not short format
    if strcmpi(format_dem_str,'i')
        try
            aps_systemcall(xyz2grd_cmd);
        catch
             fprintf('int16 does not seem to work. Let try short instead \n')
             format_dem_str='h';
             if strcmp(pixelreg,'gridline')==0
                xyz2grd_cmd = [GMT_string 'xyz2grd -R',num2str(xfirst),'/',num2str(xlast),'/',num2str(ylast),'/',num2str(yfirst),' -I',num2str(ncols),'+/',num2str(nrows),'+ ',demfile,' -Gtmp.grd -N',num2str(demnull),' -F -ZTL' format_dem_str];
             else
                xyz2grd_cmd = [GMT_string 'xyz2grd -R',num2str(xfirst),'/',num2str(xlast),'/',num2str(ylast),'/',num2str(yfirst),' -I',num2str(ncols),'+/',num2str(nrows),'+ ',demfile,' -Gtmp.grd -N',num2str(demnull),' -ZTL' format_dem_str];
             end
             aps_systemcall(xyz2grd_cmd);
        end
    else
        aps_systemcall(xyz2grd_cmd);
    end
end

% getting the information about the grid dem file

% getting a temp identified. 
% make it random such mutiple correction methods can be ran simultaneous
temp_num = round(rand(1)*10000);

y_first_new_cmd = ['echo `' GMT_string 'grdinfo tmp.grd | grep y_min`>', 'temp' num2str(temp_num)];
aps_systemcall(y_first_new_cmd);
temp = fileread(['temp' num2str(temp_num)]);
ix_y_min = findstr('y_min',temp);
ix_y_max = findstr('y_max',temp);
ix_y_end = findstr('y_inc',temp);
y_first_new = str2num(temp(ix_y_min+7:ix_y_max-2));
y_last_new = str2num(temp(ix_y_max+7:ix_y_end-2));
clear y_first_new_cmd temp ix_y_end ix_y_max ix_y_min
x_first_new_cmd = ['echo `' GMT_string 'grdinfo tmp.grd | grep x_min`>', 'temp' num2str(temp_num)];
aps_systemcall(x_first_new_cmd);
temp = fileread(['temp' num2str(temp_num)]);
ix_x_min = findstr('x_min',temp);
ix_x_max = findstr('x_max',temp);
ix_x_end = findstr('x_inc',temp);
x_first_new = str2num(temp(ix_x_min+7:ix_x_max-2));
x_last_new = str2num(temp(ix_x_max+7:ix_x_end-2));
clear x_first_new_cmd temp ix_x_end ix_x_max  ix_x_min
delete(['temp' num2str(temp_num)])

if ymin<y_first_new 
   fprintf('Your min latitude crop is outside the DEM extend, reset to the maximum \n')
   ymin = y_first_new;
end
if ymax>y_last_new
   fprintf('Your max latitude crop is outside the DEM extend, reset to the maximum \n')
   ymax = y_last_new;
end
% checking if both the region and grd are defined in the same quadrant
if xmin/abs(xmin)~=x_first_new./abs(x_first_new)
    x_first_new = x_first_new-360;
end
if xmax/abs(xmax)~=x_last_new./abs(x_last_new)
    x_last_new = x_last_new-360;
end
if xmin<x_first_new
   fprintf('Your min longitude crop is outside the DEM extend, reset to the maximum \n')
   xmin = x_first_new;
end
if xmax>x_last_new
   fprintf('Your max longitude crop is outside the DEM extend, reset to the maximum \n')
   xmax = x_last_new;
end


if exist('tmp_smp.grd','file')==2
    delete('tmp_smp.grd')
end
grdsmp_cmd = [GMT_string 'grdsample -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' -I',num2str(smpres),' tmp.grd -Gtmp_smp.grd'];
aps_systemcall(grdsmp_cmd);
grd2xyz_cmd = [GMT_string 'grd2xyz -R',num2str(xmin),'/',num2str(xmax),'/',num2str(ymin),'/',num2str(ymax),' tmp_smp.grd -bo > ',smpdem];
aps_systemcall(grd2xyz_cmd);
clear a b

%reading the dem file again for latter processing
demfid = fopen(smpdem,'r');
data_vector = fread(demfid,'double');
fclose(demfid);

% reshaping into the right n column matrix
data = reshape(data_vector,3,[])';
dem = data(:,3);
clear data data_vector


%% load the resampled DEM
if strcmp(gmt5_above,'y')
    nncols_cmd = ['echo `' GMT_string 'grdinfo tmp_smp.grd | grep n_columns | awk ''{print $NF}''`>', path_dem ,filesep ,'temp3'];
    nnrows_cmd = ['echo `' GMT_string 'grdinfo tmp_smp.grd | grep n_rows | awk ''{print $NF}''`>>', path_dem ,filesep ,'temp3'];
else
    nncols_cmd = ['echo `' GMT_string 'grdinfo tmp_smp.grd | grep nx | awk ''{print $NF}''`>', path_dem ,filesep ,'temp3'];
    nnrows_cmd = ['echo `' GMT_string 'grdinfo tmp_smp.grd | grep ny | awk ''{print $NF}''`>>', path_dem ,filesep ,'temp3'];
end
aps_systemcall(nncols_cmd);
aps_systemcall(nnrows_cmd);

DEM_info = load([path_dem ,filesep 'temp3']);
aps_systemcall(['rm ' path_dem ,filesep 'temp3']);
nncols = DEM_info(1);
nnrows = DEM_info(2);
clear DEM_info
dem =reshape(dem,nncols,nnrows)';

if fig_test ==1
    figure('name','DEM debug test');
    imagesc([xmin xmax],[ymax ymin],dem)   
    colorbar
    view(0,90)
    axis equal
    axis tight
    axis xy

    % check if this is correct
    str='';
    while ~strcmpi(str,'y') && ~strcmpi(str,'n')
        fprintf(['Does the DEM look reasonable? \n'])
        str = input('Continue? [y: for yes, n: no] \n','s');
    end
    if strcmpi(str,'n')
        error('Check the dem input file.')
    end
end
