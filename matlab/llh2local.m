function xy=llh2local(llh,origin)
%llh2local     xy=llh2local(llh,origin)
%
%Converts from longitude and latitude to local coorindates
%given an origin.  llh (lon; lat; height) and origin should
%be in decimal degrees. Note that heights are ignored and
%that xy is in km.
%
%     Copyright (C) 2015
%     Email: eedpsb@leeds.ac.uk or davidbekaert.com
%     With permission by Peter Cervelli, Jessica Murray 
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
%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Sept 7, 2000  Peter Cervelli		Original Code
%   Oct 20, 2000  Jessica Murray        Changed name from DM_llh2local to 
%                                       llh2local for use with non-DM functions;
%                                       Added to help message to clarify order
%                                       of 'llh' (i.e., lon, lat, height).
%   Dec. 6, 2000  Jessica Murray        Clarified help to show that llh 
%                                       is a column vector
%
%
%-------------------------------------------------------------

%Set ellipsoid constants (WGS84)

   a=6378137.0;
   e=0.08209443794970;

%Convert to radians
%   llh=llh*pi/180;
%   origin=origin*pi/180;
   llh=double(llh)*pi/180;
   origin=double(origin)*pi/180;

%Do the projection

   z=llh(2,:)~=0;

   dlambda=llh(1,z)-origin(1);

   M=a*((1-e^2/4-3*e^4/64-5*e^6/256)*llh(2,z) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*llh(2,z)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*llh(2,z)) - ...
        (35*e^6/3072)*sin(6*llh(2,z)));

   M0=a*((1-e^2/4-3*e^4/64-5*e^6/256)*origin(2) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*origin(2)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*origin(2)) - ...
        (35*e^6/3072)*sin(6*origin(2)));
   
   N=a./sqrt(1-e^2*sin(llh(2,z)).^2);
   E=dlambda.*sin(llh(2,z));

   xy(1,z)=N.*cot(llh(2,z)).*sin(E);
   xy(2,z)=M-M0+N.*cot(llh(2,z)).*(1-cos(E));

%Handle special case of latitude = 0

   xy(1,~z)=a*dlambda(~z);
   xy(2,~z)=-M0;

%Convert to km
   
   xy=xy/1000;
