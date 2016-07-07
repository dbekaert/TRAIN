function [x,y,z]=grdread2(file);
%GRDREAD2  Load a GMT grdfile (netcdf format)
%
% Uses NetCDF libraries to load a GMT grid file.
% Duplicates (some) functionality of the program grdread (which requires
% compilation as a mexfile-based function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% Z=GRDREAD2('filename.grd') will return the data as a matrix in Z
%
% [X,Y,Z]=GRDREAD2('filename.grd') will also return X and Y vectors
% suitable for use in Matlab commands such as IMAGE or CONTOUR.
% e.g., imagesc(X,Y,Z); axis xy
%
% Although both gridline and pixel registered grids can be read,
% pixel registration will be converted to gridline registration
% for the x- and y-vectors.
%
% See also GRDWRITE2, GRDINFO2

% CAUTION: This program currently does little error checking and makes
% some assumptions about the content and structure of NetCDF files that
% may not always be valid.  It is tested with COARDS-compliant NetCDF
% grdfiles, the standard format in GMT 4 and later, as well as GMT v3
% NetCDF formats.  It will not work with any binary grid file formats.
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
% Time-stamp: <Wed Jan  6 16:37:45 EST 2010>

% Version 1.1.1, 6-Jan-2010
% released with minor changes in documentation along with grdwrite2 and grdinfo2
% Version 1.1, 3-Dec-2009
% support for GMT v3 grids added
% Version 1.0, 29-Oct-2009
% first posted on MATLAB Central

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
  x=netcdf.getVar(ncid,0)';   
  y=netcdf.getVar(ncid,1)';
  z=netcdf.getVar(ncid,2)';
else
  if (nvars==6),                        % old (v3) GMT netCDF grid file
    [dimname, dimlen] = netcdf.inqDim(ncid,1);
    if (dimname=='xysize'),             % make sure it really is v3 netCDF
      xrange=netcdf.getVar(ncid,0)';
      yrange=netcdf.getVar(ncid,1)';
      z=netcdf.getVar(ncid,5);
      dim=netcdf.getVar(ncid,4)';
      pixel=netcdf.getAtt(ncid,5,'node_offset');
      if pixel,                         % pixel node registered
        dx=diff(xrange)/double(dim(1)); % convert int to double for division
        dy=diff(yrange)/double(dim(2));
        x=xrange(1)+dx/2:dx:xrange(2)-dx/2; % convert to gridline registered
        y=yrange(1)+dy/2:dy:yrange(2)-dy/2;
      else                              % gridline registered
        dx=diff(xrange)/double(dim(1)-1); % convert int to double for division
        dy=diff(yrange)/double(dim(2)-1);
        x=xrange(1):dx:xrange(2);
        y=yrange(1):dy:yrange(2);
      end
      z=flipud(reshape(z,dim(1),dim(2))');
    else
      error('Apparently not a GMT netCDF grid');
    end
  else
    error('Wrong number of variables in netCDF file!');
  end
end

netcdf.close(ncid)

switch nargout
  case 1,double
    varargout{1}=z;  
  case 3,
    varargout{1}=x;
    varargout{2}=y;
    varargout{3}=z;
  otherwise
    error('grdread2: Incorrect # of output arguments!');
end
