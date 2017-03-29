% function loads the APS corrections and applies it to the specified data file.
% This function assumes that you used the TRAIN2InsarApp fucntion first.
%
%     Copyright (C) 2016  Bekaert David 
%     davidbekaert.com
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
% Bekaert David 
% modifications


%%%% IFG FILE specific
filename_2cor = 'filt_topophase.unw.geo'        % data filename to be corrected. 
                                                % In case of 2-chanel data the second channel will be corrected.

%%%% APS FILE specific
% which correction needs to be converted
aps_corr = 'era';               % e.g. 'era' for ERA-I
output_prefix = '';             % if [] will take the dates, to have no prefix put to ''.


%%%% IFG STRUCTURE specific
% IFG folders are "date1-date2" or "date1_date2". 
% How are date1 and date2 defined?
date_format = 'yyyymmdd';           % this needs to be recognized by matlab.
                                    % mmm for jan format, MM is min two digit, 
                                    % while hh is hour 2 digit. You can
                                    % also include - or / if needed. e.g.
                                    % 'yyyy-mmm-dd'
                                    
                                    
%% no changes required below
% check if this is being ran from the right directory.
if exist('parms_aps.mat','file')~=2
    error('This does not look to be the right directory where you ran TRAIN')
end
% check if the file to be corrected is specified
if isempty(filename_2cor)
    error('Need to specify filename_2cor')
end

% How to add new TRAIN correction methods? Simple define below
% aps_str{}: is part of the filename_aps. Note a .geo will be added automatically if needed
if strcmpi(aps_corr,'era')
    aps_str{1} = 'era.aps';                 
    aps_str{2} = 'era_hydro.aps';
    aps_str{3} = 'era_wet.aps';
else
    error('Currently this method is not hard-coded, only few lines needed to add other method')
end
   
%% loading the csupport information and check for consistency in ".geo" extension
load('isce2train.mat','geo_str')
[temp1,temp2,temp3] =fileparts(filename_2cor);
if strcmpi(geo_str,'.geo') && ~strcmpi(temp3,'.geo')
    fprintf(['Looks like you ran train in .geo mode\nYour file to be corrected does not have a ".geo" extension\n'])
    % Check if action is needed
    repeat=1;
    while repeat==1
        action_flag= input('Do you want to update it to .geo? [y/n] ','s');
        if strcmpi(action_flag,'y')
            repeat=0;
        elseif strcmpi(action_flag,'n')
            repeat=0;
        end
    end
    if strcmpi(action_flag,'y')
        filename_2cor = [filename_2cor '.geo'];
    else
        repeat=1;
        while repeat==1
            action_flag= input('Do you want to continue? [y/n] ','s');
            if strcmpi(action_flag,'y')
                repeat=0;
            elseif strcmpi(action_flag,'n')
                repeat=0;
            end
        end
        if strcmpi(action_flag,'n')
           error('Abord by user. Update the filename_2cor') 
        end
    end
end
if strcmpi(temp3,'.geo')
    save_name = [temp2];
else
    save_name = [temp2 temp3];
end


%% writing out each APS correction in the original folder
% loading the information on the the dates as this determines where the
% data is saved
ps = load('ps2.mat');
separator = '_';
IFG_folder = [repmat(['..' filesep '..' filesep],ps.n_ifg,1)  datestr(ps.ifgday(:,1),date_format) repmat(separator,ps.n_ifg,1) datestr(ps.ifgday(:,2),date_format)];
% check if the seperation is correct
if exist(IFG_folder(1,:),'dir')~=7
    separator ='-';
    IFG_folder(IFG_folder=='_')=separator;
end
if exist(IFG_folder(1,:),'dir')~=7
    error(['Cannot recognize your IFG folders ' date_format '[_/-]' date_format])
end


% loop over the different correction methods components
for k=1:length(aps)
    fprintf(['*** ' aps{k} ':\n'])

    % loop over the different IFGS
    for k_ifg =1:ps.n_ifg
       
        % file name of the APS
        if isempty(output_prefix) && ~isstr(output_prefix)
            filename_aps = [IFG_folder(k_ifg,:) filesep  datestr(ps.ifgday(k_ifg,1),date_format) separator datestr(ps.ifgday(k_ifg,2),date_format) '_' ];
        else
            if isstr(output_prefix) && isempty(output_prefix)
                 filename_aps = [IFG_folder(k_ifg,:) filesep output_prefix(2:end) ];
            else
                filename_aps = [IFG_folder(k_ifg,:) filesep output_prefix  '_' ];
            end
        end
        filename_aps = [filename_aps  aps_str{k} geo_str];
        % the filename of the data to be corrected
        filename_2cor_data = [IFG_folder(k_ifg,:) filesep filename_2cor];
        % filename of the corrected product
        filename_out = [IFG_folder(k_ifg,:) filesep  save_name '_' aps_str{k} geo_str];

        % loading the APS data
        APS = load_isce(filename_aps);

        % Loading the data to be corrected
        [DATA] = load_isce(filename_2cor_data);
        % check the number of bands of the data, and when needed only keep second band
        n_bands = size(DATA,3);
        if n_bands>2
            DATA(:,:,3:end)=[];
        end
        if n_bands>1
            DATA(:,:,1)=[];
        end
        
        % correcting the data
        DATA = DATA-APS;
        
        % saving the data again
        DATA = reshape(DATA,[],1);
        fid = fopen(filename_out,'w');
        fwrite(fid,DATA,'float32');
        fclose(fid);
        
        % also copy the xml file from the template
        copyfile([filename_aps '.xml' ],[filename_out '.xml'])
        % update the xml file with the correct information
        set_parm_xml([filename_out '.xml'],'data_type','FLOAT');
        set_parm_xml([filename_out '.xml'],'extra_file_name','N/A');
        set_parm_xml([filename_out '.xml'],'description',['["[''corrected with APS from TRAIN using: ' aps_str{k} '.'']"]']);
        [temp1, temp2, temp3] = fileparts(filename_out);
        set_parm_xml([filename_out '.xml'],'file_name',[temp2 temp3]);
        fprintf([num2str(k_ifg) '/' num2str(ps.n_ifg) '\n'])
        clear temp1 temp2 temp3
        
        
    end
    
end