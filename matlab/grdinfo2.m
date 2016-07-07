function d = grdinfo2(file);
%GRDINFO2  Print information about a GMT grdfile (netCDF format, GMT v3 or v4)
%
% Uses NetCDF libraries to display information about a GMT grid file.
% Duplicates (some) functionality of the program grdinfo (which requires
% compilation as a mexcdf function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% GRDINFO2('file.grd') will display information about the GMT grid
% file 'file.grd' in a format similar to the gmt command grdinfo.
%
% D = GRDINFO('file.grd') will in addition return a vector containing
% (xmin, xmax, ymin, ymax, zmin, zmax, format, xinc, yinc). Format is
% 1 for pixel registration and 0 for grid node registration.
%
% See also GRDREAD2, GRDWRITE2

% This program is expected to work on any GMT netCDF format file,
% but it does not duplicate all the functionality of GMT's I/O
% library, so it will not work on all files supported by GMT.
% In particular, it will fail on binary format grdfiles.
% It is the responsibility of the user to determine whether this
% program is appropriate for any given task.
%
% For more information on GMT grid file formats, see:
% http://www.soest.hawaii.edu/gmt/gmt/doc/gmt/html/GMT_Docs/node70.html
% Details on Matlab's native netCDF capabilities are at:
% http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html

% GMT (Generic Mapping Tools, <http://gmt.soest.hawaii.edu>)
% was developed by Paul Wessel and Walter H. F. Smith

% Kelsey Jordahl
% Marymount Manhattan College
% http://marymount.mmm.edu/faculty/kjordahl/software.html

% Time-stamp: <Wed Jan  6 16:26:46 EST 2010>
 
% Version 1.1.1, 6-Jan-2010
% first released on MATLAB Central
% modification
% DB 	05/2014 	Suppress command window output

if nargin < 1,
  help(mfilename);
  return,
end

% check for appropriate Matlab version (>=7.7)
V=regexp(version,'[ \.]','split');
if (str2num(V{1})<7) | (str2num(V{1})==7 & str2num(V{2})<7),
  ver
  error('grdread2: Requires Matlab R2008b or later!');
end

ncid = netcdf.open(file, 'NC_NOWRITE');
if isempty(ncid),
  return,
end

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

if (nvars==3),                        % new (v4) GMT netCDF grid file
%  for attnum=0:(ngatts-1),
  title=getatt_clean(ncid,'title');
  conv=getatt_clean(ncid,'Conventions');
  pixel=getatt_clean(ncid,'node_offset');
  if isempty(pixel), pixel=0; end
  desc=getatt_clean(ncid,'description');
  command=getatt_clean(ncid,'history');
  vers = getatt_clean(ncid,'GMT_version');
  xname = netcdf.getAtt(ncid,0,'long_name');
  yname = netcdf.getAtt(ncid,1,'long_name');
  zname = netcdf.getAtt(ncid,2,'long_name');
  xrange = netcdf.getAtt(ncid,0,'actual_range');
  yrange = netcdf.getAtt(ncid,1,'actual_range');
  zrange = netcdf.getAtt(ncid,2,'actual_range');
  [dimname,nx]=netcdf.inqDim(ncid,0);
  [dimname,ny]=netcdf.inqDim(ncid,1);
else
  if (nvars==6),                        % old (v3) GMT netCDF grid file
    [dimname, dimlen] = netcdf.inqDim(ncid,1);
    if (dimname=='xysize'),             % make sure it really is v3 netCDF
      title=getatt_clean(ncid,'title');
      command=getatt_clean(ncid,'source');
      conv=[];
      desc=[]; zrange=[0 0]; vers='3.x format file';
      xrange=netcdf.getVar(ncid,0)';
      xname=netcdf.getAtt(ncid,0,'units');
      yrange=netcdf.getVar(ncid,1)';
      yname=netcdf.getAtt(ncid,1,'units');
      zrange=netcdf.getVar(ncid,2)';
      zname=netcdf.getAtt(ncid,2,'units');
%      z=netcdf.getVar(ncid,5);
      dim=netcdf.getVar(ncid,4)';
      nx=dim(1); ny=dim(2);
      zname='z';
      pixel=netcdf.getAtt(ncid,5,'node_offset');
    else
      error('Apparently not a GMT netCDF grid');
    end
  else
    error('Wrong number of variables in netCDF file!');
  end
end

if pixel,                         % pixel node registered
  dx=diff(xrange)/double(nx); % convert int to double for division
  dy=diff(yrange)/double(ny);
else                              % gridline registered
  dx=diff(xrange)/double(nx-1); % convert int to double for division
  dy=diff(yrange)/double(ny-1);
end

%disp(['Title: ' title]);
%disp(['Conventions: ' conv]);
%disp(['GMT version: ' vers]);
%disp(['Command: ' command]);
%disp(['Remark: ' desc]);
%if pixel,
%  disp('Pixel node registration used');
%else
%  disp('Gridline node registration used');
%end
%disp(['x_min: ' num2str(xrange(1)) ' x_max: ' num2str(xrange(2)) ...
%                   ' name: ' xname])
%disp(['y_min: ' num2str(yrange(1)) ' y_max: ' num2str(yrange(2)) ...
%                   ' name: ' yname])
%disp(['z_min: ' num2str(zrange(1)) ' z_max: ' num2str(zrange(2)) ...
%                   ' name: ' zname])

switch nargout
  case 1,
   % need to make sure that zrange and pixel are cast into double precision
   d=[xrange(1) xrange(2) yrange(1) yrange(2) double(zrange(1)) double(zrange(2)) double(pixel) dx dy];
end

function val = getatt_clean(ncid,attname)
%
% Call netcdf.getAtt when not sure whether global attribute 'attname' exists.
% Trap error and return empty vector in that case.
%
eval(['val = netcdf.getAtt(ncid,netcdf.getConstant(''NC_GLOBAL''),''' attname ''');'],'val=[];');
