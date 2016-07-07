function [] = aps_powerlaw_update_band(band_number)
% function that updates ther variable in tca2_sb or tca2 for the selected
% powerlaw band.
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
% By Bekaert David -- University of Leeds
%
% modifications
stamps_processed = getparm_aps('stamps_processed');

if strcmp(stamps_processed,'y')
    load psver
else
    psver = 2; 
end

save_path = [pwd filesep]; 

% file names of the output data
apsname = [save_path filesep 'tca' num2str(psver) '.mat'];
apssbname = [save_path filesep 'tca_sb' num2str(psver) '.mat'];
apsbandsname = [save_path filesep 'tca_bands' num2str(psver) '.mat'];
apsbandssbname = [save_path filesep 'tca_bands_sb' num2str(psver) '.mat'];



% loading the data variable
if strcmp(stamps_processed,'y')
    if strcmp(getparm('small_baseline_flag'),'y')
       aps_file = apssbname;
       aps_band_file = apsbandssbname;
    else
        aps_file = apsname;
        aps_band_file = apsbandsname;
    end
else
    aps_file = apsname;
    aps_band_file = apsbandsname;
end


% loading the band data
load(aps_band_file,'K_tropo_powerlaw_bands','K_tropo_powerlaw_or','ph_tropo_powerlaw_bands','ph_tropo_powerlaw_or');
if band_number==0
    K_tropo_powerlaw = K_tropo_powerlaw_or;
    ph_tropo_powerlaw = ph_tropo_powerlaw_or;
    fprintf('Updated the powerlaw results back to the orginal (all bands combined)\n')
else
    K_tropo_powerlaw = K_tropo_powerlaw_bands(:,:,band_number);
    ph_tropo_powerlaw = ph_tropo_powerlaw_bands(:,:,band_number);
    sptial_bands = getparm_aps('powerlaw_spatial_bands');
    fprintf(['Updated the powerlaw to band ' num2str(round(sptial_bands(band_number,1)/1000*10)/10) '-' num2str(round(sptial_bands(band_number,2)/1000*10)/10) ' km' '\n'])
end

% saving the data angain
if exist(aps_file,'file')==2
    save(aps_file,'-append','ph_tropo_powerlaw','K_tropo_powerlaw')
else
    save(aps_file,'ph_tropo_powerlaw','K_tropo_powerlaw')    
end

