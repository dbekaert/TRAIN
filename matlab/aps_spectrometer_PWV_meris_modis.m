function [calibration_matrix,corr_vector,fraction_vector,dates] = aps_spectrometer_PWV_meris_modis
% aps_spectrometer_PWV_meris_modis(batchfile)
% Scipt to load meris and modis PWV data and compare each SAR data for the
% cloud free pixels.
%
% Bekaert David -- University of Leeds
% 
% modifications:
% DB    08/2014     Initial codings
% DB    03/2015     Add bootstrapping to get the uncertancy of the linear
%                   mapping between PWV of MERIS and MODIS.


save_dir= 'aps_m_M';
n_boot_runs = 1200;

curdir = pwd;

fontsize = 15;

% loading the date information  
meris_datapath = getparm_aps('meris_datapath',1);
modis_datapath = getparm_aps('modis_datapath',1);
stamps_processed = getparm_aps('stamps_processed',1);
if strcmp(stamps_processed,'y')
   ps = load(getparm_aps('ll_matfile',1));
   dates = ps.day;
   fprintf('Stamps processed structure \n')
else
    ifgday_matfile = getparm_aps('ifgday_matfile',1);
    dates = load(ifgday_matfile);
    dates = dates.ifgday;
    dates = reshape(dates,[],1);
    dates = unique(dates);
end



%start loop here to calculate atmos correction for each date
ndates = length(dates);
fprintf('Starting the computation for each SAR date \n')
calibration_matrix = NaN([ndates 2]);            % this is the linear correction matrix [a b] between MERIS = a*MODIS+b (PWV);
corr_vector = NaN([ndates 1]);
fraction_vector =NaN([ndates 1]);
std_vector=NaN([ndates 2]);
calibration_matrix_boot = NaN([n_boot_runs 2 ndates]);

if exist([curdir filesep 'figures'],'dir')~=7
    mkdir([curdir filesep 'figures']);
end


for n = 1:ndates
    date_str = datestr(dates(n,1),'yyyymmdd');

    meris_filename_PWV_nointerp = [meris_datapath filesep date_str filesep date_str '_ZPWV_nointerp.xyz'];
    modis_filename_PWV_nointerp = [modis_datapath filesep date_str filesep date_str '_ZPWV_nointerp.xyz'];

    % checking if there is actual meris data for this date, if not just
    % leave NaN's in the matrix.
    if exist(meris_filename_PWV_nointerp,'file') ==2 && exist(modis_filename_PWV_nointerp,'file') ==2
      
        [xyz_meris_PWV,temp] = load_meris_SAR(meris_filename_PWV_nointerp);
        [xyz_modis_PWV,temp] = load_meris_SAR(modis_filename_PWV_nointerp);
        clear temp

        % keep only those points in common
        lonlat_modis = xyz_modis_PWV(:,1:2);
        lonlat_meris = xyz_meris_PWV(:,1:2);
        
        % number of original grid points, this also includes water!
        n_original = size(lonlat_meris,1);
        
        [lonlat, ix_modis, ix_meris ]= intersect(lonlat_modis,lonlat_meris,'rows','stable');
        
        PWV_modis=xyz_modis_PWV(ix_modis,3);
        PWV_meris=xyz_meris_PWV(ix_meris,3);
        
        ix_nan = sum(isnan([PWV_modis PWV_meris]),2)>0;
        lonlat(ix_nan,:)=[];
        PWV_modis(ix_nan,:)=[];
        PWV_meris(ix_nan,:)=[];
        
        if length(PWV_meris)>100
            % number of gridpoitns left
            n_used = size(lonlat,1);
            fraction_vector(n,1) = n_used./n_original;
            
            % getting the starts
            A = [PWV_modis  ones(size(PWV_meris)) ];
            calibration_matrix(n,:) = lscov(A,PWV_meris);
            corr_vector(n,1) = corr(PWV_modis,PWV_meris);
            
            
            % getting the uncertaincy through bootstrapping - Can just use
            % the powerlaw bootstrap function as its not technique dependent
            [std_vector(n,:),calibration_matrix_boot(:,:,n)]= aps_powerlaw_bootstrap(A,PWV_meris,n_boot_runs);
            
            
            limits = [floor(min([PWV_modis;PWV_meris]))   ceil(max([PWV_modis;PWV_meris]))];
            h1 = figure('name',['MODIS/MERIS PWV comparison ' date_str]);
            plot(limits, limits,'k--')
            hold on
            plot(limits,calibration_matrix(n,1).*limits+calibration_matrix(n,2),'r-','linewidth',2)
            hold on
            plot(PWV_modis,PWV_meris,'k.')
            xlim(limits)
            ylim(limits)
            hold on
            plot(limits,calibration_matrix(n,1).*limits+calibration_matrix(n,2),'r-','linewidth',2)
            hold on
            plot(limits, limits,'k--')
            xlabel('MODIS PWV','fontsize',fontsize)
            ylabel('MERIS PWV','fontsize',fontsize)
            title({[date_str ':'],['PWV_{MERIS} = ' num2str(round(calibration_matrix(n,1)*100)/100) ' [+-' num2str(round(std_vector(n,1)*100)/100)  '] *PWV_{MODIS} + ' num2str(round(calibration_matrix(n,2)*100)/100) ' [+-' num2str(round(std_vector(n,2)*100)/100) ']' ]},'fontsize',fontsize)
            legend('1-to-1 relation','Linear fit',0)
            set(gca,'fontsize',fontsize)         
            print(h1,'-depsc',[curdir filesep 'figures' filesep 'PWV_MERIS_MODIS_' date_str '.eps'])            
            print(h1,'-dpng',[curdir filesep 'figures' filesep 'PWV_MERIS_MODIS_' date_str '.png'])
            close(h1)
        end
        fprintf([num2str(n) ' completed out of ' num2str(ndates) '\n']) 
    else
        fprintf([date_str ': no comparison\n']) 
    end
    
    
end

save([modis_datapath filesep 'MODIS_calibration.mat'],'calibration_matrix','corr_vector','fraction_vector','dates','std_vector','calibration_matrix_boot');
