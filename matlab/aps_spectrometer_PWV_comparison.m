function [] = aps_spectrometer_PWV_comparison(start_step,end_step)
% [] = aps_spectrometer_PWV_comparison(start_step,end_step)
% Compares the Precipitable Water Vapor estimates from MERIS and MODIS.
% MODIS is know to over-estimate the PWV. This code allows to get an
% indication of this calibration factor for MODIS.
%
% By David Bekaert - University of Leeds
% August 2014
%
% Modifications:
% DB	08/2014 	Convert the program to a function - Initial codings



%% loading the StaMPS variables
curdir = pwd;
modis_datapath = getparm_aps('modis_datapath',1);


if start_step==1 
    % Getting PWV estiamted from both MERIS and MODIS
    fprintf('\n\nStep 1: Extract MERIS/MODIS Percipitable Water Vapor for individual dates\n\n')
        
    %% doing MERIS    
    fprintf('Doing MERIS...\n')
    % generating a file list of the files to be processed:
    meris_datapath = getparm_aps('meris_datapath',1);
    command_str = 'echo files > meris_batch_file.txt';
    command_str2 = ['ls -d ' meris_datapath filesep '2*' filesep 'MER_RR_*reprojected.tif >> meris_batch_file.txt'];

    [temp, temp1] = system(command_str);
    [temp, temp1] = system(command_str2);
    clear temp temp1
        
    % running the MERIS extraction on the file list
    aps_spectrometer_PWV_meris('meris_batch_file.txt')
    
    %% doing MODIS
    fprintf('\nDoing MODIS...\n')
    % generating a file list of the files to be processed:
    command_str = 'echo files > modis_batch_file.txt';
    command_str2 = ['ls -d ' modis_datapath filesep '2*' filesep 'OSCAR_Modis_*.grd >> modis_batch_file.txt'];

    [temp, temp1] = system(command_str);
    [temp, temp1] = system(command_str2);
    clear temp temp1
                
    % running the MERIS computation on the file list
    aps_spectrometer_PWV_modis('modis_batch_file.txt')
    
end

if start_step<=2 && end_step >=2
    fprintf('\n\nStep 2: Comparing MERIS/MODIS PWV for SAR dates\n\n')
    [calibration_matrix,corr_vector,fraction_vector,dates] = aps_spectrometer_PWV_meris_modis;   
end
if start_step<=3 && end_step >=3
    fprintf('\n\nStep 3: Updating MODIS PWV for SAR dates\n\n')
    curdir = pwd;
    load([modis_datapath filesep 'MODIS_calibration.mat']);
    
    
    calibration_matrix_new = calibration_matrix;
    corr_vector_new = corr_vector;
    fraction_vector_new = fraction_vector;
    
    % getting those were no comparison was made for
    ix = isnan(corr_vector);
    
    % checking if the correlation is ok. This indicates point spread around
    % linear relationship.
    ix1 = corr_vector<0.5;
    if sum(ix1)>0
        fprintf([num2str(sum(ix1)) ' was/were was rejected because there is not a clear linear relationship \n']) 
        calibration_matrix_new(ix1,:) = NaN;
        corr_vector_new(ix1,:) = NaN;
        fraction_vector_new(ix1,:) = NaN;
    end
    
    % checking if there is too few data but it to 10 % of all data
    ix2 = fraction_vector_new<0.1;
    if sum(ix2)>0
        fprintf([num2str(sum(ix2)) ' was/were rejected because did not meet the 10perc data threshold \n']) 
        calibration_matrix_new(ix2,:) = NaN;
        corr_vector_new(ix2,:) = NaN;
        fraction_vector_new(ix2,:) = NaN;
    end
    ix1 = ix1+ix2;
    clear ix2
    
    calibration_info = [corr_vector corr_vector_new];
    if sum(isnan(corr_vector_new))==length(corr_vector_new)
        fprintf('No SAR dates for which the calibration factor could be estimated \n')
    else
        fprintf([num2str(sum(~isnan(corr_vector_new))) ' out of ' num2str(length(corr_vector_new)) ' have a calibration estimated \n']);
        fprintf(['Substituting the mean for the missing dates\n']);
        mean_calibration = nanmean(calibration_matrix_new,1);
        calibration_matrix_new(isnan(corr_vector_new),:) = repmat(nanmean(calibration_matrix_new,1),sum(isnan(corr_vector_new)),1);
        calibration_matrix = calibration_matrix_new;
        
        for counter=1:length(corr_vector_new)
            if counter==1
                fprintf(['\nmean correction: MERIS PWV = ' num2str(round(mean_calibration(1)*100)/100) ' * MODIS PWV + ' num2str(round(mean_calibration(2)*100)/100) '\n'])
                fprintf('----------------------------------------------------------------------\n')
                fprintf('%s\t %s\t %s\t %s\n','Date','MERIS PWV = cor * MODIS PWV', ' + offset','info')
            end
            if ix(counter)>0
                fprintf('%s\t %3.2f\t\t  %3.2f\t\t %s\n',datestr(dates(counter),'yyyymmdd'),calibration_matrix(counter,1),calibration_matrix(counter,2),' no comparison, fixed with mean')                
            elseif ix1(counter)>0
                fprintf('%s\t %3.2f\t\t  %3.2f\t\t %s\n',datestr(dates(counter),'yyyymmdd'),calibration_matrix(counter,1),calibration_matrix(counter,2),' fixed with mean')                
            else
                fprintf('%s\t %3.2f\t\t  %3.2f\n',datestr(dates(counter),'yyyymmdd'),calibration_matrix(counter,1),calibration_matrix(counter,2))
            end
            if counter==length(corr_vector_new)
                fprintf('----------------------------------------------------------------------\n\n')
            end
        end
        

        % check with user if the files need to be updated with the new
        % calibration factor
        str_recal = '';
        while strcmpi(str_recal,'y')~=1 && strcmpi(str_recal,'n')~=1
            str_recal = input(['Do you want to use this to re-calibrate MODIS? [y/n] \n'],'s');
        end
        
        % no reprocessing, do you want to visualize it?
        if strcmpi(str_recal,'y')
            % generating a file list of the files to be processed:
            modis_datapath = getparm_aps('modis_datapath');
            command_str = 'echo files > modis_batch_file.txt';
            command_str2 = ['ls -d ' modis_datapath filesep '2*' filesep 'OSCAR_Modis_*.grd >> modis_batch_file.txt'];
            [temp, temp1] = system(command_str);
            [temp, temp1] = system(command_str2);
            clear temp temp1

                     
            % check if correction is using mean or on individual dates
            str_recal = '';
            while strcmpi(str_recal,'y')~=1 && strcmpi(str_recal,'n')~=1
                str_recal = input(['Use the average calibration factor (if n, uses a different calibration factor for each SAR date)? [y/n] \n'],'s');
            end
            
            % use average
            if strcmpi(str_recal,'y')
                calibration_matrix = repmat(mean_calibration,size(calibration_matrix,1),1);
            end
            

            % running the MODIS computation on the file list
            aps_modis_recalibrate('modis_batch_file.txt',calibration_matrix,dates,calibration_info)
            fprintf('To use this data put modis_recalibrated to ''y'' using setparm_aps, and run aps_modis. \n')
        end
    end
end


cd(curdir)



