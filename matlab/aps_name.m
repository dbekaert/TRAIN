function [technique_str, save_str] = aps_name(aps_str)
% fucntion that return the string for a legend based on the input
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

if strcmp(aps_str,'a_p')
    technique_str = 'Power-law (w+h)';
    save_str = 'Powerlaw_WH';
elseif strcmp(aps_str,'a_m')
    technique_str = 'Meris (w)';
    save_str = 'Meris_W';
elseif strcmp(aps_str,'a_mi')
    technique_str = 'Meris (w)';
    save_str = 'Meris_W';
elseif strcmp(aps_str,'a_m+a_eh')
    technique_str = 'Meris (w) + ECMWF (h)';
    save_str = 'Meris_W_ECMWF_H';
elseif strcmp(aps_str,'a_m+a_wh')
    technique_str = 'Meris (w) + WRF (h)';
    save_str = 'Meris_W_WRF_H';
elseif strcmp(aps_str,'a_mi+a_eh')
    technique_str = 'Meris (w) + ECMWF (h)';
    save_str = 'Meris_W_ECMWF_H';
elseif strcmp(aps_str,'a_mi+a_wh')
    technique_str = 'Meris (w) + WRF (h)';
    save_str = 'Meris_W_WRF_H';
elseif strcmp(aps_str,'a_M')
    technique_str = 'Modis (w)';
    save_str = 'Modis_W';
elseif strcmp(aps_str,'a_MI')
    technique_str = 'Modis (w)';  
    save_str = 'Modis_W';
elseif strcmp(aps_str,'a_M+a_eh')
    technique_str = 'Modis (w) + ECMWF (h)';
    save_str = 'Modis_W_ECMWF_H';
elseif strcmp(aps_str,'a_M+a_wh')
    technique_str = 'Modis (w) + WRF (h)';
    save_str = 'Modis_W_WRF_H';
elseif strcmp(aps_str,'a_MI+a_eh')
    technique_str = 'Modis (w) + ECMWF (h)';
    save_str = 'Modis_W_ECMWF_H';
elseif strcmp(aps_str,'a_MI+a_wh')
    technique_str = 'Modis (w) + WRF (h)';   
    save_str = 'Modis_W_WRF_H';    
elseif strcmp(aps_str,'a_e')
    technique_str = 'ECMWF (w+h)';  
    save_str = 'ECMWF_WH';    
elseif strcmp(aps_str,'a_eh')
    technique_str = 'ECMWF (h)';  
    save_str = 'ECMWF_H';    
elseif strcmp(aps_str,'a_ew')
    technique_str = 'ECMWF (w)';
    save_str = 'ECMWF_W';    
elseif strcmp(aps_str,'a_w')
    technique_str = 'WRF (w+h)';   
    save_str = 'WRF_WH';    
elseif strcmp(aps_str,'a_wh')
    technique_str = 'WRF (h)';   
    save_str = 'WRF_H';    
elseif strcmp(aps_str,'a_ww')
    technique_str = 'WRF (w)'; 
    save_str = 'WRF_W';    
elseif strcmp(aps_str,'a_l')
    technique_str = 'Linear (w+h)';    
    save_str = 'Linear_WH';    
elseif strcmp(aps_str,'a_lman')
    technique_str = 'Linear (w+h)';   
    save_str = 'LinearMAN_WH'; 
elseif strcmp(aps_str,'a_RM')
    technique_str = 'Modis recal (w)';
    save_str = 'ModisRecal_W';
elseif strcmp(aps_str,'a_RMI')
    technique_str = 'Modis recal (w)';  
    save_str = 'ModisRecal_W';
elseif strcmp(aps_str,'a_RM+a_eh')
    technique_str = 'Modis recal (w) + ECMWF (h)';
    save_str = 'ModisRecal_W_ECMWF_H';
elseif strcmp(aps_str,'a_RMI+a_eh')
    technique_str = 'Modis recal (w) + ECMWF (h)';
    save_str = 'ModisRecal_W_ECMWF_H';
else
    aps_str
    error('Do not know this type')
end
        

