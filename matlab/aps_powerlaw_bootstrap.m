function [coeff_std,coeff_vector]= aps_powerlaw_bootstrap(A,y,bood_num)
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
% By Bekaert David 

% iterations
for k=1:bood_num

    % bootstrapping positions
    ix = ceil(rand(size(y,1),1)*size(y,1));       % gives position of the earlier points to be used of i1
    A_new = A;
    y_new = y;
    A_new = A_new(ix,:);                      % gives the position wrt all data points
    y_new = y_new(ix,:);
  
   
    coeff = lscov(A_new,y_new);
    % correct for the introduced scaling
    coeff_vector(k,:) = coeff;
end    
coeff_std = std(coeff_vector,0);