function [] = aps_merra_files(orderflag,merra_model)
%     Copyright (C) 2016  Bekaert David
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
% By David Bekaert - April 2016
%
%
% MERRA data overview at website:
% http://disc.sci.gsfc.nasa.gov/uui/datasets/GES_DISC_MAI6NPANA_V5.2.0/summary?keywords=%22MERRA%22
% MERRA has data from 1979 till 
% download data page is at: ftp://goldsmr3.sci.gsfc.nasa.gov/data/s4pa/MERRA/MAI6NPANA.5.2.0/
% download data page for subset is different:
% If the year is < 1993:
% http://goldsmr3.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAI6NPANA.5.2.0%2F/YEAR/%2F/MONTH/%2FMERRA100.prod.assim.inst6_3d_ana_Np./DATE/.hdf&FORMAT=SERGLw&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T/HOUR/%3A00%3A00%2F1979-01-01T/HOUR/%3A00%3A00&LABEL=MERRA100.prod.assim.inst6_3d_ana_Np./DATE/.SUB.hdf&FLAGS=&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Ch%2Ct%2Cqv
% If the year is < 2001:
% http://goldsmr3.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAI6NPANA.5.2.0%2F/YEAR/%2F/MONTH/%2FMERRA200.prod.assim.inst6_3d_ana_Np./DATE/.hdf&FORMAT=SERGLw&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T/HOUR/%3A00%3A00%2F1979-01-01T/HOUR/%3A00%3A00&LABEL=MERRA200.prod.assim.inst6_3d_ana_Np./DATE/.SUB.hdf&FLAGS=&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Ch%2Ct%2Cqv
% If the year is >=2001:
% http://goldsmr3.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAI6NPANA.5.2.0%2F/YEAR/%2F/MONTH/%2FMERRA300.prod.assim.inst6_3d_ana_Np./DATE/.hdf&FORMAT=SERGLw&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T/HOUR/%3A00%3A00%2F1979-01-01T/HOUR/%3A00%3A00&LABEL=MERRA300.prod.assim.inst6_3d_ana_Np./DATE/.SUB.hdf&FLAGS=&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Ch%2Ct%2Cqv

% MERRA2 data overview at website:
% http://disc.sci.gsfc.nasa.gov/uui/datasets/GES_DISC_M2I6NPANA_V5.12.4/summary?keywords=%22MERRA%22
% MERRA2 is continuation of MERRA and has data from 1979 till current
% downlaod data page is at: ftp://goldsmr5.sci.gsfc.nasa.gov/data/s4pa/MERRA2/M2I6NPANA.5.12.4/
% modifications:
% SSS    10/2016 MERRA credentials now passed
% SSS    10/2016 Url change for MERRA-1 data processed for Summer 2010
% DB     11/2017 Give the end date for the MERRA model
% DB     11/2017 Attempt to fix URL for merra download
% DB     02/2018 Add more modular approach between models

if nargin<1 || isempty(orderflag)
    orderflag=0;
end
if nargin<2
    merra_model='merra';
end

%% Loading dataset specific paramters
% read username /n password for MERRA credentials, which are to be provided in '~/.merrapass'
if orderflag==1
    fileID = fopen('~/.merrapass','r');
    if fileID==-1
       error('~/.merrapass containing credentials does not exist or is not readable, also add NASA GESDISC DATA ARCHIVE in your applications to access the archive')
    end
    permis= textscan(fileID,'%s');
    fclose(fileID);
    usern=permis{1}{1};
    pass=permis{1}{2};
    clear fileID permis;
end

% getting the variables from the parms_aps file
workdir = pwd;

% datapath specific
merra_datapath = getparm_aps('merra_datapath',1);
if isempty(merra_datapath)
    error('please specify merra_datapath')
end

% loading the data
stamps_processed = getparm_aps('stamps_processed');
UTC_sat =  getparm_aps('UTC_sat');
ifgday_matfile = getparm_aps('ifgday_matfile');
ifgs_dates = load(ifgday_matfile);
ifg_dates = ifgs_dates.ifgday;
if ~strcmp(stamps_processed,'y')
    ifg_dates = datenum(num2str(ifg_dates),'yyyymmdd');
end


%% Compute based on Satellite pass which weather model outputs that will be used
[time_before,time_after, date_before, date_after,f_before,f_after] = aps_weather_model_times(merra_model,ifg_dates,UTC_sat);
time_vector = [time_before ; time_after];
date_vector = [date_before ;  date_after];
clear f_before f_after

%% downloading the data
% number of files that needs to be downloaded
n_files = size(date_vector,1);

% extend the BBOX a bit larger than the data bbox
[S,N,W,E] = aps_weather_model_crop;

% now generate the full paths where the data can be downloaded
wget_string = repmat('wget ',n_files,1);
if strcmpi(merra_model,'merra')

    % year specific url
    datasetnumber = 200.*ones([n_files 1]);
    datasetnumber(str2num(date_vector(:,1:4))<1993)=100;
    datasetnumber(str2num(date_vector(:,1:4))>=2001)=300;
    datasetnumber(str2num(date_vector(:,1:8))>=20100601 & str2num(date_vector(:,1:8))<=20100831)=301; %%%SSS 10/16: Url change for MERRA-1 data processed for Summer 2010.
    datasetnumber_str = num2str(datasetnumber);
    
    % the URL of the file              
    filesTOdownlaod = [repmat('http://goldsmr3.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA%2FMAI6NPANA.5.2.0%2F',n_files,1) date_vector(:,1:4) repmat('%2F',n_files,1) date_vector(:,5:6) repmat('%2FMERRA',n_files,1) datasetnumber_str repmat('.prod.assim.inst6_3d_ana_Np.',n_files,1) date_vector(:,1:8) repmat('.hdf&FORMAT=bmM0Lw&BBOX=',n_files,1) repmat(num2str(S),n_files,1)  repmat('%2C',n_files,1) repmat(num2str(W),n_files,1) repmat('%2C',n_files,1) repmat(num2str(N),n_files,1) repmat('%2C',n_files,1) repmat(num2str(E),n_files,1) repmat('&TIME=1979-01-01T',n_files,1) time_vector(:,1:2) repmat('%3A00%3A00%2F1979-01-01T',n_files,1) time_vector(:,1:2) repmat('%3A00%3A00&LABEL=MERRA300.prod.assim.inst6_3d_ana_Np.',n_files,1) date_vector(:,1:8) repmat('.SUB.nc4&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_MERRA&VERSION=1.02&LAYERS=&VARIABLES=ps%2Ch%2Ct%2Cqv',n_files,1)];

    
    % the filename of the downloaded file as to be stored
    downloadFILEname = [repmat([merra_datapath filesep],n_files,1) date_vector(:,1:8) repmat(filesep,n_files,1)  repmat('MERRA_',n_files,1) date_vector(:,1:8)  repmat('_',n_files,1) time_vector(:,1:2)  repmat('.hdf',n_files,1)];
    downloadFILEname = [repmat([merra_datapath filesep],n_files,1) date_vector(:,1:8) repmat(filesep,n_files,1)  repmat('MERRA_',n_files,1) date_vector(:,1:8)  repmat('_',n_files,1) time_vector(:,1:2)  repmat('.nc4',n_files,1)];
       
    % give warning about end of the model
    ix_drop = datenum(date_vector,'yyyymmdd')>datenum('20160229','yyyymmdd');
    if sum(ix_drop)>1
        fprintf('Note MERRA model only has outputs till 20160229 \n')
        fprintf([ num2str(sum(ix_drop)) '/' num2str(length(ix_drop)) ' dates do not have MERRA outputs \n'])
    end
    downloadFILEname(ix_drop,:)=[];
    filesTOdownlaod(ix_drop,:)=[];
    
elseif strcmpi(merra_model,'merra2')   
   
    % year specific url
    datasetnumber = 200.*ones([n_files 1]);
    datasetnumber(str2num(date_vector(:,1:4))<1992)=100; %%%SSS 7/2016: Should be <1992 instead of 1993.
    datasetnumber(str2num(date_vector(:,1:4))>=2001)=300;
    datasetnumber(str2num(date_vector(:,1:4))>=2011)=400;
    datasetnumber_str = num2str(datasetnumber);
    
    % the URL of the file
    filesTOdownlaod = [repmat('http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA2%2FM2I6NPANA.5.12.4%2F',n_files,1) date_vector(:,1:4) repmat('%2F',n_files,1) date_vector(:,5:6) repmat('%2FMERRA2_',n_files,1) datasetnumber_str repmat('.inst6_3d_ana_Np.',n_files,1) date_vector(:,1:8) repmat('.nc4&FORMAT=bmM0Yy8&BBOX=',n_files,1) repmat(num2str(S),n_files,1)  repmat('%2C',n_files,1) repmat(num2str(W),n_files,1) repmat('%2C',n_files,1) repmat(num2str(N),n_files,1) repmat('%2C',n_files,1) repmat(num2str(E),n_files,1) repmat('&TIME=1979-01-01T',n_files,1) time_vector(:,1:2) repmat('%3A00%3A00%2F1979-01-01T',n_files,1) time_vector(:,1:2) repmat('%3A00%3A00&LABEL=svc_MERRA2_',n_files,1) datasetnumber_str repmat('.inst6_3d_ana_Np.',n_files,1) date_vector(:,1:8) repmat('.nc4&FLAGS=&SHORTNAME=M2I6NPANA&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES=',n_files,1)];
    
    % the filename of the downloaded file as to be stored
    downloadFILEname = [repmat([merra_datapath filesep],n_files,1) date_vector(:,1:8) repmat(filesep,n_files,1)  repmat('MERRA2_',n_files,1) date_vector(:,1:8)  repmat('_',n_files,1) time_vector(:,1:2)  repmat('.nc4',n_files,1)];

else
    error('Does not exist. Either MERRA or MERRA2 model')
end


% order the data and download it.
overwrite_flag=-1;
% update number of files for MERRA 1 could have dropped few to download
n_files = size(downloadFILEname,1);
if orderflag==1
    % not that some symbols from the url needs escape char to get through using wget.
    for k=1:n_files
        if exist(downloadFILEname(k,:),'file')==0
            [path,downloadfile,file_ext] = fileparts(downloadFILEname(k,:));
            if exist(path,'dir')~=7
                mkdir(path);
            end
            fprintf(['Downloading: ' downloadfile file_ext '\n'])
            try
                pause(5); 
                pass_to_cmd=['wget --user ',usern,' --password ', pass,' ''',filesTOdownlaod(k,:),'''',' -O ',downloadFILEname(k,:)];
                [a,b] = system(pass_to_cmd);
                clear a b pass_to_cmd;
            catch ME
                fprintf('File not found \n')
                fprintf('Matlab failure message: %s\n', ME.message);
                clear ME;
            end
        else
            [path,downloadfile,file_ext] = fileparts(downloadFILEname(k,:));
            if overwrite_flag==-1
                str='';
                while ~strcmpi(str,'y') && ~strcmpi(str,'n') 
                    fprintf(['Do you want to overwrite existing files? \n'])  
                    str = input('[y: for yes, n: no] \n','s');
                end
                if strcmpi(str,'n')
                    overwrite_flag=0;
                else
                    overwrite_flag=1;
                end
            end
            % check if the files need to be overwritten
            if overwrite_flag==1
                delete(downloadFILEname(k,:))
                fprintf(['Downloading: ' downloadfile file_ext '\n'])
                try
                    pause(5); 
                    pass_to_cmd=['wget --user ',usern,' --password ', pass,' ''',filesTOdownlaod(k,:),'''',' -O ',downloadFILEname(k,:)];
                    [a,b] = system(pass_to_cmd);
                    clear a b pass_to_cmd;
                catch ME
                    fprintf('File not found \n')
                    fprintf('Matlab failure message: %s\n', ME.message);
                    clear ME;
                end
            elseif overwrite_flag==0
                fprintf('File %s has already been downloaded \n',[downloadfile file_ext ])
            end
        end
    end
else    
    % outputing this information to a file
    fid = fopen([upper(merra_model) '_files.txt'],'w');
    fprintf(['Required ' upper(merra_model) ' files for all interferograms \n'])
    for k=1:size(filesTOdownlaod,1)
        fprintf('%s \n',[filesTOdownlaod(k,:)]);
        fprintf(fid,'%s \n',[filesTOdownlaod(k,:)]);
    end
    fclose(fid);
end
