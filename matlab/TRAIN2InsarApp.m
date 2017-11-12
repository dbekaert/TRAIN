% function that re-writes data back into radar or geo coordinates for ISCE
% this function assumes that you used the InsarApp2TRAIN function.
% IF you did not use this function originally, note that there are some assumption 
% on how the data was originally stored. E.g. when you put a radar coordinate
% data into TRAIN it was assumed to be done according to data_train = reshape(data,[],1);
% Also there is information read from a support file 'isce2train.mat'.
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
% 8/10/2016     DB      Update the xml information to reflect TRIAN
%                       parameters and compatible with mdx.py

% which correction needs to be converted
aps_corr = 'era';               % e.g. 'era' for ERA-I
output_prefix = '';             % if [] will take the dates, to have no prefix put to ''.

% IFG folders are "date1-date2" or "date1_date2". 
% How are date1 and date2 defined?
date_format = 'yyyymmdd';           % this needs to be recognized by matlab.
                                    % mmm for jan format, MM is min two digit, 
                                    % while hh is hour 2 digit. You can
                                    % also include - or / if needed. e.g.
                                    % 'yyyy-mmm-dd'
                                    
                                    
%% no changes required below
% check if this is being ran from the right directory.
if exist('tca_sb2.mat','file')~=2 && exist('tca2.mat','file')~=2
    if exist('parms_aps.mat','file')~=2
        error('You have not ran any TRAIN correction yet')
    else
        error('This does not look to be the right directory where you ran TRAIN')
    end
end
if strcmpi(getparm_aps('small_baseline_flag'),'y')
   tca = 'tca_sb2.mat';
else
    tca = 'tca2.mat';           % this should normally not occur when using the InsarApp2TRAIN script.
end


% How to add new TRAIN correction methods? Simple define below
% aps_str{}: is part of the final filename. Note a .geo will be added automatically if needed
% aps{}: the variable name as know in tca2.mat or tca_sb2.mat for TRAIN.
if strcmpi(aps_corr,'era')
    aps{1} = 'ph_tropo_era';
    aps{2} = 'ph_tropo_era_hydro';
    aps{3} = 'ph_tropo_era_wet';
    aps_str{1} = 'era.aps';                 
    aps_str{2} = 'era_hydro.aps';
    aps_str{3} = 'era_wet.aps';
else
    error('Currently this method is not hard-coded, only few lines needed to add other method')
end


    
%% loading the data specific information
load('isce2train.mat','MASKED_flag','geo_str','templatexmlfile')
% if masked then load the IND directly 
if strcmpi(MASKED_flag,'y')
    load('isce2train.mat','IND','WIDTH','LENGTH')
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
    eval(['data = load(''' tca ''',''' aps{k} ''');'])
    try 
        eval(['data =  data.' aps{k} ';'])
    catch
       error(['You have not calculated the ' aps_corr ' correction yet...']) 
    end
    fprintf(['*** ' aps{k} ':\n'])

    % loop over the different IFGS
    for k_ifg =1:ps.n_ifg
        % converting from vector back to radar coordinates
        if strcmpi(MASKED_flag,'y')
            [vector_data] = Dvec2Dmat(data(:,k_ifg),WIDTH,LENGTH,IND,0);
            vector_data = reshape(vector_data,[],1);
            % fprintf('There is a more efficient way to do the masking\n')
        else
            vector_data = data(:,k_ifg);
        end

        % save file names
        if isempty(output_prefix) && ~isstr(output_prefix)
            filename = [IFG_folder(k_ifg,:) filesep  datestr(ps.ifgday(k_ifg,1),date_format) separator datestr(ps.ifgday(k_ifg,2),date_format) '_' ];
        else
            if isstr(output_prefix) && isempty(output_prefix)
                 filename = [IFG_folder(k_ifg,:) filesep output_prefix(2:end) ];
            else
                filename = [IFG_folder(k_ifg,:) filesep output_prefix  '_' ];
            end
        end
        filename = [filename  aps_str{k} geo_str];
        fid = fopen(filename,'w');
        fwrite(fid,vector_data,'float32');
        fclose(fid);
        % also copy the xml file from the template
        copyfile(['..' filesep '..' filesep templatexmlfile],[filename '.xml'])
        
        % update the xml file with the correct information
        set_parm_xml([filename '.xml'],'data_type','FLOAT');
        set_parm_xml([filename '.xml'],'extra_file_name','N/A');
        set_parm_xml([filename '.xml'],'description',['["[''Estimated APS from TRAIN using: ' aps_str{k} '.'']"]']);
        [temp1, temp2, temp3] = fileparts(filename);
        set_parm_xml([filename '.xml'],'file_name',[temp2 temp3]);
        fprintf([num2str(k_ifg) '/' num2str(ps.n_ifg) '\n'])
        clear temp1 temp2 temp3
        
    end
end

fprintf(['\n\n****\nIn case you clean the TRAIN processing folder, make sure to keep:\n - "ps2.mat" \n - "parms_aps.mat" \n - "isce2train.mat" \nThese are used by ApplyCor2InsarApp.m function\n*****\n\n'])