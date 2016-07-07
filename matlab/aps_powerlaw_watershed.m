function [mountain_ridge] = aps_powerlaw_watershed()
% [] = aps_powerlaw_watershed()
% Script to compute a watershed of the DEM. This can be used to define
% patches in such a way the patches edges are placed at the ridges.
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
%
% Bekaert David - Leeds University 2014
%
% modifications
% 07/2014       DB  Save information in the tca support variable


fontsize = 15;
resolution = 800;       % put dem resolution to 100 m


% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed');

% getting the variables from the parms_aps file
stamps_processed = getparm_aps('stamps_processed');
ll_matfile = getparm_aps('ll_matfile');
lonlatfile= load(ll_matfile);
xlims = [min(lonlatfile.lonlat(:,1)) max(lonlatfile.lonlat(:,1))];
ylims =  [min(lonlatfile.lonlat(:,2)) max(lonlatfile.lonlat(:,2))];
ix_convexhull = convhull(lonlatfile.lonlat(:,1),lonlatfile.lonlat(:,2));


%% Compute and resample DEM
[dem,xmin,xmax,ymin,ymax,smpres,nncols,nnrows] = get_DEM;

% smoothing the DEM
dem_smooth = medfilt2(dem, [6 6]);
dem_smooth(dem_smooth==0)=max(max(dem_smooth));


h1=figure;
imagesc([xmin xmax],[ymax ymin],dem_smooth)   
cc=colorbar
view(0,90)
axis equal
axis tight   
colormap(flipud(gray))
axis xy
xlabel(cc,'[m]','fontsize',fontsize)
ylabel(cc,'Topography','fontsize',fontsize)
set(gca,'fontsize',fontsize)

deminfo.xmin = xmin;
deminfo.xmax = xmax;
deminfo.ymax = ymax;
deminfo.ymin = ymin;
deminfo.dem_smooth = dem_smooth;


hold on
plot(lonlatfile.lonlat(ix_convexhull,1),lonlatfile.lonlat(ix_convexhull,2),'g-','linewidth',2)


counter=1;
continue_flag=1;
while continue_flag

    if counter==1
      fprintf('\n\nStart selecting your points, press enter once complete \n')
      fprintf('You will be asked for a soft or a hard ridge.\nWindows constrain eachother across soft ridges, while they don''t over hard ridges \n')
    end

  mountain_ridge{counter}=ginput;


  % ask if this is a hard ridge or not
  str='';
  while ~strcmpi(str,'y') && ~strcmpi(str,'n') 
    str = input('Is this a hard ridge? (y) yes, or (n) no \n','s');
  end
  if strcmpi(str,'y')
      hard_ridge_flag(counter)=1;
  else
      hard_ridge_flag(counter)=0;
  end

  if hard_ridge_flag(counter)==0;
      % plotting the line on the figure
      plot(mountain_ridge{counter}(:,1),mountain_ridge{counter}(:,2),'r--','linewidth',2)
  else
      % plotting the line on the figure
      plot(mountain_ridge{counter}(:,1),mountain_ridge{counter}(:,2),'r-','linewidth',2)
  end
  
  str='';
  while ~strcmpi(str,'n') && ~strcmpi(str,'q') 
    str = input('n: for next line clasification, q: to stop  \n','s');
  end
  if strcmpi(str,'q')
      continue_flag=0;
      InSAR_convexhull = lonlatfile.lonlat(ix_convexhull,:);
      powerlaw_ridges.mountain_ridge = mountain_ridge;
      powerlaw_ridges.InSAR_convexhull = InSAR_convexhull;
      powerlaw_ridges.deminfo = deminfo;
      powerlaw_ridges.hard_ridge_flag = hard_ridge_flag;
          
          
      if exist('tca_support.mat','file')==2
          save('tca_support.mat','-append','powerlaw_ridges','InSAR_convexhull')
      else
          save('tca_support.mat','powerlaw_ridges','InSAR_convexhull')        
      end
      print(h1,'-dpng',['aps_p' filesep 'mountain_ridges.png'])
      print(h1,'-depsc',['aps_p' filesep 'mountain_ridges.eps'])
  else
      counter = counter+1;
  end      
end



