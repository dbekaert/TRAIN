function [value,parmname]=getparm_aps(parmname,printflag)
%GETPARM_APS get parameter value from parms_lwts.mat
% GETPARM_APS(PARMNAME) 
% Only enough characters of PARMNAME to make it unique need be typed
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
% Based on script by Andy Hooper (StaMPS)
% Modified for the aps toolbox by David Bekaert - University of Leeds - 2013
% This script allows for non-StaMPS structured processed data.

parms_default_aps

if nargin<2
    printflag=0;
end

parmfile='parms_aps.mat';

if exist(['.' filesep parmfile],'file')
    parms=load(parmfile);
elseif exist(['..' filesep parmfile],'file')
    parmfile=['..' filesep parmfile];
    parms = load(parmfile);
else
    error([parmfile ' not found'])
end


if nargin < 1
    disp(orderfields(parms))
else
    parmnum=strmatch(parmname,fieldnames(parms)); 
    if length(parmnum)>1
        error(['Parameter ',parmname,'* is not unique'])
    elseif isempty(parmnum)
        parmname=[];
        value=[];
    else
        parmnames=fieldnames(parms);
        parmname=parmnames{parmnum};
        value=getfield(parms,parmname);
    end
    if printflag~=0
        if isnumeric(value)
            fprintf(['   PARM: %s=',repmat('%g ',1,200)],parmname,value)
            fprintf('\n')
        else
            fprintf('   PARM: %s=''%s''\n',parmname,value)
        end
    end
end



    
    

