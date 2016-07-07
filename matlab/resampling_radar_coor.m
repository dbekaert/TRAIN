function [ij,xy,phase_matrix]=resampling_radar_coor(multi_look_range,phase_in,ps)
% Multi-look inputted data based on the radar coordiantes and and the
% multi-look factor for the range direction. The azimuth direction follows
% based on the aspect ratio. As output the multi-looked data is given
% together with a the new positions XY in a local frame and the radar
% coordinates IJ.
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
% Bekaert David



if nargin<1
    multi_look_range = 20;
end

ar = ps.ar;
ps.ij_new = ps.ij;
ps.ij_new(:,3) = ceil(ps.ij(:,3)./max(ps.ij(:,3))*multi_look_range);
ps.ij_new(:,2) = ceil(ps.ij(:,2)./ar./max(ps.ij(:,3))*multi_look_range);


ij = NaN([1 3]);
xy = NaN([1 3]);
lonlat = NaN([1 2]);
phase = NaN([1 size(phase_in,2)]);


counter=1;
n_grid = max(ps.ij_new(:,2))*max(ps.ij_new(:,3));
h = waitbar(0,'Resampling data, Please wait...');
for i=1:max(ps.ij_new(:,2))
    
    for k=1:max(ps.ij_new(:,3))
        ix = find(ps.ij_new(:,2)==i & ps.ij_new(:,3)==k);
        if isempty(ix)~=1            
            ij(counter,1) = counter;
            xy(counter,1) = counter;
            ij(counter,2:3) = [i k];
            xy(counter,2:3) = mean(ps.xy(ix,2:3),1); 
            lonlat(counter,1:2) = mean(ps.lonlat(ix,:),1);


            A = ones([size(ix,1) 1]);
            phase(counter,:) = lscov(A,phase_in(ix,:) ); 


            counter = counter +1;
            waitbar(counter/n_grid)
        end
        
    end
end
close(h)

% reconstruct the matrix in radar coordiantes
phase_matrix =  NaN([max(ij(:,2:3))  size(phase_in,2) ]);

for k=1:ij(end,1)
    phase_matrix(ij(k,2),ij(k,3),:)=phase(k,:);
end


bperp = ps.bperp;
day = ps.day;
day_ix = ps.day_ix;
master_day = ps.master_day;
master_ix = ps.master_ix;
n_ifg = ps.n_ifg;
n_image = ps.n_image;
n_ps = size(lonlat,1);









