function [] = aps_era_ECMWF_Python(timing,wheatherregstr)
% script for the creation of download python script files for ERA-I model
% downloads monthly data in splitted times
% input str: timing in the format 'YYYYMMDDHH'
%
%     Copyright (C) 2015  
%     With permission of Hannes Bathke
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
% By Hannes Bathke -- January 2014
%
% modifications:
% DB    02/2014 Integrate in the APS toolbox
% DB    03/2014 Auto download included

year = timing(1:4);
month = timing(5:6);
day = timing(7:8);
htime = timing(9:10);

%mfield = calendar(str2double(year),str2double(month));
%nday = max(max(mfield));
datestr = [year,'-',month,'-',day]; 

pythonsc_path = [timing,'.py'];
fid = fopen(pythonsc_path,'w');
fprintf(fid,'#!/usr/bin/env python\n');
fprintf(fid,'\n');
fprintf(fid,'from ecmwfapi import ECMWFDataServer\n');
fprintf(fid,'\n');
fprintf(fid,'# To run this example, you need an API key\n');
fprintf(fid,'# available from https://api.ecmwf.int/v1/key/\n');
fprintf(fid,'\n');
fprintf(fid,'server = ECMWFDataServer()\n');
fprintf(fid,'\n');
fprintf(fid,'server.retrieve({\n');
fprintf(fid,'          ''dataset'' : "interim",\n');
fprintf(fid,'          ''class''   : "ei",\n');
fprintf(fid,'          ''type''    : "an",\n');
fprintf(fid,'          ''stream''  : "oper",\n');
fprintf(fid,'          ''levtype'' : "pl",\n');
fprintf(fid,'          ''levelist'': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",\n');
fprintf(fid,'          ''param''   : "129.128/130.128/157.128/246.128",\n');
fprintf(fid,'          ''date''    : "%s",\n',datestr);
fprintf(fid,'          ''time''    : "%s",\n',htime);
fprintf(fid,'          ''step''    : "0",\n');
fprintf(fid,'          ''format''  : "netcdf",\n');
fprintf(fid,'          ''resol''   : "auto",\n');
fprintf(fid,'          ''area''    : "%s",\n',wheatherregstr);
fprintf(fid,'          ''grid''    : "0.75/0.75",\n');
fprintf(fid,['          ''target''  : "ggap' num2str(timing(1:10)) '00.nc",\n']);
fprintf(fid,'                })\n');


