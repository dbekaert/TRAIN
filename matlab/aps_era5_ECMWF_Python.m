function [] = aps_era5_ECMWF_Python(timing,wheatherregstr)
% script for the creation of download python script files for ERA5 model
% downloads data for given date and time
% input str: timing in the format 'YYYYMMDDHH'
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
% By Bekaert Davide -- June 2017
%
% modifications:

year = timing(1:4);
month = timing(5:6);
day = timing(7:8);
htime = timing(9:10);
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
fprintf(fid,'          ''dataset'' : "era5_test",\n');
fprintf(fid,'          ''type''    : "an",\n');
fprintf(fid,'          ''stream''  : "oper",\n');
fprintf(fid,'          ''levtype'' : "pl",\n');
fprintf(fid,'          ''param''   : "129.128/130.128/157.128/246.128",\n');
fprintf(fid,'          ''date''    : "%s",\n',datestr);
fprintf(fid,'          ''time''    : "%s",\n',htime);
fprintf(fid,'          ''format''  : "netcdf",\n');
fprintf(fid,'          ''grid''    : "0.25/0.25",\n');


fprintf(fid,'          ''step''    : "0",\n');
fprintf(fid,'          ''levelist'': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",\n');
fprintf(fid,'          ''resol''   : "auto",\n');
fprintf(fid,'          ''area''    : "%s",\n',wheatherregstr);


fprintf(fid,['          ''target''  : "ggap' num2str(timing(1:10)) '00.nc",\n']);
fprintf(fid,'                })\n');


