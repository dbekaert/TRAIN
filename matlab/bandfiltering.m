function [] = bandfiltering(z_grid,x_res,y_res,spatial_bands,save_path,ifg_based_correction,n_degree_butterworth,norm_filter_flag)
% function that computes the band pass filtered data usigng 2D FFT and a
% butterwurth function. Make sure the input data is on a regular grid and
% has no NaN values. There is no output passed in this function. Instead
% the bandfiltered data is saved for each dataset individually, e.g.
% 'dataset_1.mat' for the first dataset and increasing for the others. 
% input:
% z_grid                    The data specified as a grid, i.e. a matrix.
%                           Additonal datasets can be specified by 
%                           increasing the third dimention.
% x_res                     The resolution in x-direction in m of the grid
% y_res                     The resolution in y-direction in m of the grid
% spatial_bands             The spatial band that need to be filtered in m, 
%                           specified as a 2 column matrix [lower upper].
%                           Multiple band filters can be specified by 
%                           increasing the number or rows. 
% Optional inputs:
% save_path                 The path were the output data will be saved. default
%                           is the current directory
% n_degree_butterworth      The degree of the butterwurth filter, by
%                           default this is set to be 3.
% norm_filter_flag          Normalisation of the butterwurth filter is done
%                           by default. This is to scope with the issue
%                           when the extremes of the bandfitler are too
%                           close to eachother causing only a partial
%                           amplitude passing in the selected band.
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
% By Bekaert David - December 2012
% modifications:
% DB	02/2013		Include the save_path variable
% DB    03/2013     Include warnings for edge effects and to large spatial
%                   bandfilters.
% DB    03/2013     Include 1D filtering to cope with limitation of narrow 
%                   datasets  
% DB    05/2013     Include option to crop out a region

% optional inputs for checking results:
% figure properties
fontsize = 15;
plot_figures = 0;          % plot the figures for the n_plot_dataset dataset
save_fig =0;                % when 1 save figures in the figures folder in tha aps_p folder.
n_plot_dataset=1;           % plot the figures for this band dataset. 
                            % Topography is first then the interferograms
check_flag = 0;             % some additional figures being generated when turn on
mirror_flag = 1;            % when 1, mirror the dataset such filtering effects
                            % are reduced. The mirror is based on the
                            % largest spatial filter wavelength.
                            
max_perc_mirror_2D = 100;    % filtersize with respect to the dataset spatial 
                            % dimension in percentage till what 2D
                            % filtering is allowed. When larger 1D
                            % filtering is performed in the larger
                            % dimension of X or Y. 

                            
warning_flag_perc = 10;     % Output a warning when more no mirroring is done,
                            % and the percentage of half the filter length 
                            % with respect to the maximum dimension is mirrored.
                            
warning_mirror_flag_perc = 50;  
                            % Output a warning when more than this percentage 
                            % with respect to the maximum dimension is mirrored.
                            % Mirroring is done by half the maximum filter size.

                            
% getting the data from the parms_aps file
crop_flag = getparm_aps('crop_flag');
 % setting the function defaults
if nargin<6
    ifg_based_correction='n';
end
if nargin<7
    n_degree_butterworth = [];
end
if nargin<8
    norm_filter_flag = [];
end
if isempty(n_degree_butterworth)==1
   n_degree_butterworth=3; 
end
if isempty(norm_filter_flag)==1
   norm_filter_flag=1; 
end
if nargin<5 || isempty(save_path)==1
   save_path = './';
end
fprintf(['***Bandfiltering***\n'])
fprintf(['Using butterworth filter degree: ' num2str(n_degree_butterworth) '\n'])
if norm_filter_flag==1
    fprintf(['Normalise the butterworth filter. \n'])
end

norm_filter_flag=0
n_degree_butterworth = 10                                            
% number of datasets
n_datasets = size(z_grid,3);
n_band_filters = size(spatial_bands,1);
if strcmp(ifg_based_correction,'y')
    h_ifg_number = n_datasets/2;
else
    h_ifg_number=1;
end                            
                            
if size(spatial_bands,1)==1 && size(spatial_bands,2)==2 && spatial_bands(1,1)==0 && spatial_bands(1,2)==inf
    % this is a band filter of the whole image
        for k=1:n_datasets

            if k<=h_ifg_number && h_ifg_number~=1
                % these are interferograms
                save_name = ['bandfilter_regular_hgt_ifg_' num2str(k) '.mat'];
            elseif k<=h_ifg_number && h_ifg_number==1
                % these are the heights
                save_name = 'bandfilter_regular_hgt.mat';
            else
                % thse are interferograms
                 save_name = ['bandfilter_regular_ifg_' num2str(k-h_ifg_number) '.mat'];
            end
            data_band_out(:,:) = z_grid(:,:,k);
            clear data_band
            
            dimension_filter=NaN;
            % saving the bandfiltered data for eacht dataset seperately
            save([save_path, filesep, save_name],'data_band_out','x_res','y_res','spatial_bands','dimension_filter')
            clear data_band_out 
        end
else

    % Starting the bandfiltering code
    if plot_figures==1   
        if n_plot_dataset<=h_ifg_number && h_ifg_number~=1
            % thse are interferograms
            save_folder_str=['aps_p' filesep 'fig_bandfilter_hgt_ifg' num2str(n_plot_dataset)];
        elseif n_plot_dataset<=h_ifg_number && h_ifg_number==1
            % these are the heights
            save_folder_str=['aps_p' filesep 'fig_bandfilter_hgt'];
        else
            % these are interferograms
            save_folder_str=['aps_p' filesep 'fig_bandfilter_ifg_' num2str(n_plot_dataset-h_ifg_number)];
        end

    %     if n_plot_dataset==1
    %         save_folder_str=['aps_p' filesep 'fig_bandfilter_hgt'];
    %     else
    %         save_folder_str=['aps_p' filesep 'fig_bandfilter_ifg_' num2str(n_plot_dataset-1)];
    %     end
        if exist(save_folder_str,'dir')~=7
           mkdir(save_folder_str) 
        end
    end

    % original gridsize
    data_rows_or = size(z_grid,1);
    data_columns_or = size(z_grid,2);

    % investigate in the percentage of padding that is being performed
    % below the percentage of the maximum band filter with respect to the dimension is computed
    x_samples_overlap = ceil(max(spatial_bands,[],2)./x_res/2);
    x_perc_overlap = x_samples_overlap./data_columns_or*100;
    y_samples_overlap = ceil(max(spatial_bands,[],2)./y_res/2);
    y_perc_overlap = y_samples_overlap./data_rows_or*100;

    % output a warning in case the padding is to big or when the the filter
    % effect hit the warning percentage as set at the start of the code.
    if mirror_flag==1
       warning_perc = warning_mirror_flag_perc ;
    else
       warning_perc = warning_flag_perc ;    
    end
    ix_x = find(x_perc_overlap>warning_perc);
    ix_y = find(y_perc_overlap>warning_perc);
    ix = intersect(ix_x,ix_y);
    clear ix_x ix_y
    if isempty(ix)~=1 
        new_linestr = repmat('\n',length(ix),1);        
        outputstr = [num2str(spatial_bands(ix,:)) new_linestr];
        if mirror_flag==1
            fprintf(['***Warning: The following band filters are mirrored about ',num2str(warning_perc),' perc of the maximum dimension: \n'])
            for kk=1:size(outputstr,1)
                fprintf(outputstr(kk,:)) ;
            end
            clear outputstr kk
            fprintf(['This migth introduce arctifacts. \n'])
            fprintf(['It is recommended to limit to smaller spatial bandwidths. \n'])

        else
            fprintf(['***Warning: The following band filters are above ',num2str(warning_perc),' perc of the maximum dimension: \n'])
            for kk=1:size(outputstr,1)
                fprintf(outputstr(kk,:)) ;
            end
            clear outputstr kk
            fprintf(['Likely edge effects are introduced. \n'])
            fprintf(['Turn the mirror flag on and/or limit to smaller spatial bandwidths. \n\n'])
        end
    end
    clear ix

    % output information on which bandfilters are replaced by a 1D filter to
    % reduce edge effects from a limiting dimension
    % x/y_perc_overlap represents the maximum bandfilter size with respect to
    % the x/y dimension given as a percentage
    ix_x = find(x_perc_overlap>max_perc_mirror_2D);
    ix_y = find(y_perc_overlap>max_perc_mirror_2D);
    
    if isempty(ix_x)~=1 || isempty(ix_y)~=1
        x_n_cases = length(ix_x);
        y_n_cases = length(ix_y);
        if x_n_cases>y_n_cases
             fprintf(['X-direction appears limited for band filtering. \n']) 
             fprintf(['1D (Y-direction) band filtering is performed for the following bands: \n']) 
             new_linestr = repmat('\n',length(ix_x),1);
             outputstr = [num2str(spatial_bands(ix_x,:)) new_linestr];
             for kk=1:size(outputstr,1)
                 fprintf(outputstr(kk,:)) ;
             end

             % setting the variables for the 1D filtering
             filter_1D_Y = 1;
             filter_1D_X = 0;
             ix_1D_filter = ix_x;
             clear outputstr new_linestr kk ix_y ix_x

        elseif y_n_cases>x_n_cases
             fprintf(['Y-direction appears limited for band filtering. \n']) 
             fprintf(['1D (X-direction) band filtering is performed for the following bands: \n']) 
             new_linestr = repmat('\n',length(ix_y),1);
             outputstr = [num2str(spatial_bands(ix_y,:)) new_linestr];
             for kk=1:size(outputstr,1)
                 fprintf(outputstr(kk,:)) ;
             end                  

             % setting the variables for the 1D filtering
             filter_1D_Y = 0;
             filter_1D_X = 1;
             ix_1D_filter = ix_y;
             clear outputstr new_linestr kk ix_x ix_y

        else
             fprintf(['Both dataset appears as limited for band filtering. \n']) 
             fprintf(['2D band filtering is performed, but check the following bands: \n']) 
             new_linestr = repmat('\n',length(ix_y),1);
             outputstr = [num2str(spatial_bands(ix_y,:)) new_linestr];
             for kk=1:size(outputstr,1)
                 fprintf(outputstr(kk,:)) ;
             end   

             % setting the variables for the 1D filtering
             filter_1D_Y = 0;
             filter_1D_X = 0;
             ix_1D_filter = [];
             clear outputstr new_linestr kk ix_x ix_y
        end
    else
         filter_1D_Y = 0;
         filter_1D_X = 0;
         ix_1D_filter = [];
    end


    if plot_figures==1
        data_original_image = z_grid(:,:,n_plot_dataset);
    end


    % mirror the edges of the grid to reduce filter effects at the edges
    % get the maximum size of the spatial filters
    if mirror_flag == 1
        fprintf('Perform mirroring to reduce filtering effects on edges \n')
        fprintf('By half the maximum filter length at each edge. \n')

        % maximum spatial wavelength
        max_band = max(max(spatial_bands));
        if max_band == inf
           max_band = 100000;        % extend the grid to a maximum of 100 km in case one goes to infinite 
        end
        n_mirror_x = ceil(max_band./x_res*1.5);
        n_mirror_y = ceil(max_band./y_res*1.5);

        for k=1:n_datasets
            % padding the grid symmetric 
            z_grid_new(:,:,k) = padarray(z_grid(:,:,k),[n_mirror_y n_mirror_x],'symmetric');
            % plotting an intermediate figure when requested
            if plot_figures==1 && k==n_plot_dataset
                h1= figure('name','Data symmetric padded for largest filter');
                imagesc(z_grid_new(:,:,k))
                axis equal
                axis tight

                % saving of the figure when requested
                if save_fig==1
                    fig_save_name = [save_folder_str filesep 'original_data_mirrored.eps'];
                    set(h1,'PaperPositionMode','auto')
                    print(h1,'-depsc','-r150',fig_save_name)
                    clear fig_save_name
                    close(h1)
                end
                clear h1 
            end
        end
        clear z_grid
        z_grid = z_grid_new;
    else
       fprintf('No mirroring performed. Filter artifacts will become more persistent for larger wavelengths! \n') 
    end


    % size of the dataset
    data_rows = size(z_grid,1);
    data_columns = size(z_grid,2);

    % Sampling frequency follows from the resolution
    fs_rows = 1/y_res;          % rows sampling frequency [1/m]
    fs_columns = 1/x_res;   	% columns sampling frequency [1/m]

    % rows and columns of the data such they are a power of 2
    % this will speed up the fft and will autmoatically padd 
    % the data matrix with zeros where needed
    data_rows_new = 2.^nextpow2(data_rows);
    data_columns_new = 2.^nextpow2(data_columns);

    if plot_figures==1
        if mirror_flag==1
            % computation of the axis extremes based on the resolution given
            X_lims_fig = [0 data_columns_or*x_res];     % axis limits in [m]
            Y_lims_fig = [0 data_rows_or*y_res];        % axis limits in [m] 
        else
            X_lims_fig = [0 data_columns*x_res];        % axis limits in [m]
            Y_lims_fig = [0 data_rows*y_res];           % axis limits in [m] 
        end
        % computation of the axis extremes based on the resolution given 
        % Asummed to be symmetric padded. In case no padding is done, this will
        % be equal to the orginal dataset
        X_lims_fig_syix_1D_filterm = [0 data_columns*x_res];        % axis limits in [m]
        Y_lims_fig_sym = [0 data_rows*y_res];           % axis limits in [m]
    end


    % Loop over each dataset. Invert to the freq domain once and reinvert for
    % each bandfilter. i.e. There will be a loop over the different
    % bandfilters. For each dataset the data is being saved in a matfile.
    for k=1:n_datasets
        if k==n_plot_dataset && plot_figures==1
            plot_fig_flag = 1;
        else
            plot_fig_flag = 0;
        end

        fprintf(['Progress: ' num2str(k) '/' num2str(n_datasets) ,' done \n'])
        % when chosen visualise the dataset
        if plot_fig_flag==1
            h1 = figure('name','Original image');
            imagesc(X_lims_fig,Y_lims_fig,data_original_image)
            cc = colorbar;
            axis xy
            axis equal
            axis tight
            set(gca,'fontsize',fontsize)
            xlabel('Distance [m]','fontsize',fontsize)
            ylabel('Distance [m]','fontsize',fontsize)
            title('Original image','fontsize',fontsize)
            colorlimits = get(cc,'YLim');

            % saving of the figure when requested
            if save_fig==1
                fig_save_name = [save_folder_str filesep 'input.eps'];
                set(h1,'PaperPositionMode','auto')
                print(h1,'-depsc','-r150',fig_save_name)
                clear fig_save_name
                close(h1)
            end
            clear h1 


        end

        % converting to the frequency domain by FFT
        data_freq_complex = fft2(z_grid(:,:,k),data_rows_new,data_columns_new);	

        % shifting the spectrum around zero freq
        % matlab shifts with 1 pixel off for even number columns and rows
        data_freq_complex = fftshift(data_freq_complex);
        data_freq_complex_new = [data_freq_complex  data_freq_complex(:,1)];
        data_freq_complex_new = [data_freq_complex_new ; data_freq_complex_new(1,:)];
        clear data_freq_complex

        % Frequency domain: sampling resolution follows from the Nyquist requency (spampling freq/samples)
        % freq resolution of the bins in the freq domain
        fs_rows = 1/data_rows_new*1/y_res;		
        fs_columns = 1/data_columns_new*1/x_res;

        % Computing the frequencies for the axis in the freq domain figures.
        freq_rows_fig_vector = [-(data_rows_new-1)/2-0.5:1:(data_rows_new-1)/2+0.5].*fs_rows;
        freq_columns_fig_vector= [-(data_columns_new-1)/2-0.5:1:(data_columns_new-1)/2+0.5].*fs_columns;

        if plot_fig_flag==1  && check_flag==1
            figure('name','Spectrum centralised around zero freq')
            imagesc(freq_columns_fig_vector,freq_rows_fig_vector,abs(data_freq_complex_new))
            xlabel('Xfreq [1/m]','fontsize',fontsize)
            ylabel('Yfreq [1/m]','fontsize',fontsize)
            title('Spectrum centralised around zero freq','fontsize',fontsize)
            colorbar
            set(gca,'fontsize',fontsize)
            axis equal
            axis tight

            figure('name','Phase centralised around zero freq')
            imagesc(freq_columns_fig_vector,freq_rows_fig_vector,angle(data_freq_complex_new))
            xlabel('Xfreq [1/m]','fontsize',fontsize)
            ylabel('Yfreq [1/m]','fontsize',fontsize)
            title('Phase centralised around zero freq','fontsize',fontsize)
            colorbar
            set(gca,'fontsize',fontsize)
            axis equal
            axis tight
        end



        %% Computing the grid of frequencies
        freq_rows_fig_matrix = repmat(freq_rows_fig_vector',1,data_columns_new+1);
        freq_columns_fig_matrix = repmat(freq_columns_fig_vector,data_rows_new+1,1);
        % Set the frequencies for 1D or 2D filtering
        if filter_1D_X==1 
            % 1D X-direction filter
            freq_fig_matrix_1D = freq_columns_fig_matrix;
            if plot_fig_flag==1 && check_flag==1
                figure('name','Frequency grid 1D X-direction filter')
                imagesc(freq_fig_matrix_1D)
                xlabel('X','fontsize',fontsize)
                ylabel('Y','fontsize',fontsize)
                cc = colorbar;
                set(gca,'fontsize',fontsize)
                xlabel(cc,'Freq [1/m]','fontsize',fontsize)
            end
        elseif filter_1D_Y==1
            % 1D Y-direction filter
            freq_fig_matrix_1D = freq_rows_fig_matrix;
            if plot_fig_flag==1 && check_flag==1
                figure('name','Frequency grid 1D Y-direction filter')
                imagesc(freq_fig_matrix_1D)
                xlabel('X','fontsize',fontsize)
                ylabel('Y','fontsize',fontsize)
                cc = colorbar;
                set(gca,'fontsize',fontsize)
                xlabel(cc,'Freq [1/m]','fontsize',fontsize)
            end
        end
        % 2D filter
        freq_fig_matrix_2D = sqrt(freq_rows_fig_matrix.^2+freq_columns_fig_matrix.^2);
        if plot_fig_flag==1 && check_flag==1
            figure('name','Frequency grid')
            imagesc(freq_fig_matrix_2D)
            xlabel('X','fontsize',fontsize)
            ylabel('Y','fontsize',fontsize)
            cc = colorbar;
            set(gca,'fontsize',fontsize)
            xlabel(cc,'Freq [1/m]','fontsize',fontsize)
        end


        %% Loop over the different bandfilters
        % initialisation of the output data_band variable
        % This is after the symmetric padding has been removed again.
        data_band_out = NaN([data_rows_or data_columns_or n_band_filters]);
        dimension_filter = NaN([n_band_filters 1]);
        for kk=1:n_band_filters
            % output to the screen
            fprintf(['Spatial bandfilter: ',num2str(spatial_bands(kk,1)), '\t - ',num2str(spatial_bands(kk,2)),' \t m \t'])

            % convert spatial wavelength band to spatial frequency band
            f_band = 1./sort(spatial_bands(kk,:));	% frequency band [1/m]


            % Checking if the band filter is set as 1D or 2D filter
            filter_1D = find(kk==ix_1D_filter);
            if isempty(filter_1D)~=1
                % do 1D filtering
                % computing the function of the bandfilter
                H_low = 1./(1+(freq_fig_matrix_1D./f_band(1)).^(2*n_degree_butterworth));
                H_high = 1./(1+(freq_fig_matrix_1D./f_band(2)).^(2*n_degree_butterworth));
                dimension_filter(kk) = 1;
                fprintf(['(1D filter)\n'])
            else
                % do 2D filtering
                % computing the function of the bandfilter
                H_low = 1./(1+(freq_fig_matrix_2D./f_band(1)).^(2*n_degree_butterworth));
                H_high = 1./(1+(freq_fig_matrix_2D./f_band(2)).^(2*n_degree_butterworth));
                dimension_filter(kk) = 2;
                fprintf(['(2D filter)\n'])
            end

            % Correct for the case of a filter having inf in its range
            ix_center_spectrum = (size(H_high)-1)./2+1;
            if isnan(H_high(ix_center_spectrum(1),ix_center_spectrum(2)))
                H_high(ix_center_spectrum(1),ix_center_spectrum(2))=0;
            end
            if isnan(H_low(ix_center_spectrum(1),ix_center_spectrum(2)))
                H_low(ix_center_spectrum(1),ix_center_spectrum(2))=0;
            end
            % Computing the band pass filter
            H_Butterworth =  (H_low - H_high);
            H_Butterworth_norm = H_Butterworth./max(max(H_Butterworth));

            % plotting when requested
            if plot_fig_flag==1
                % plotting the filter around zero for the rows
                h1 = figure('name','Rows freq filter');
                subplot(4,1,1)
                plot(freq_rows_fig_vector,H_low(:,ix_center_spectrum(2)),'b.-')
                ylim([0 1])
                title(['freq : ' num2str(f_band(1)) ])
                subplot(4,1,2)
                plot(freq_rows_fig_vector,H_high(:,ix_center_spectrum(2)),'b.-')
                ylim([0 1]) 
                title(['freq : ' num2str(f_band(2)) ])
                subplot(4,1,3)
                plot(freq_rows_fig_vector,H_Butterworth(:,ix_center_spectrum(2)),'b.-')
                ylim([0 1])
                title(['Band freq' ])
                subplot(4,1,4)
                plot(freq_rows_fig_vector,H_Butterworth_norm(:,ix_center_spectrum(2)),'b.-')
                ylim([0 1])
                title(['Norm band freq' ])

                % saving of the figure when requested
                if save_fig==1
                    fig_save_name = [save_folder_str filesep 'bandfilter_' num2str(spatial_bands(kk,1)) '_'  num2str(spatial_bands(kk,2)) 'm.eps'];
                    set(h1,'PaperPositionMode','auto')
                    print(h1,'-depsc','-r150',fig_save_name)
                    clear fig_save_name
                    close(h1)
                end
                clear h1 

            end

            if norm_filter_flag ==1 
                H_Butterworth = H_Butterworth_norm;
            end

            % multiplying filter in frequency domain.
            % data_freq_complex_new = real(data_freq_complex_new).*H_Butterworth+1i.*angle(data_freq_complex_new);
            data_freq_complex_new_band = data_freq_complex_new.*H_Butterworth;
            if plot_fig_flag==1  && check_flag==1
                figure('name','New spectrum centralised around zero freq')
                imagesc(freq_columns_fig_vector,freq_rows_fig_vector,abs(data_freq_complex_new_band))
                xlabel('Xfreq [1/m]','fontsize',fontsize)
                ylabel('Yfreq [1/m]','fontsize',fontsize)
                title('Spectrum centralised around zero freq','fontsize',fontsize)
                colorbar
                set(gca,'fontsize',fontsize)
                axis equal
                axis tight

                figure('name','New phase centralised around zero freq')
                imagesc(freq_columns_fig_vector,freq_rows_fig_vector,angle(data_freq_complex_new_band))
                xlabel('Xfreq [1/m]','fontsize',fontsize)
                ylabel('Yfreq [1/m]','fontsize',fontsize)
                title('Phase centralised around zero freq','fontsize',fontsize)
                colorbar
                set(gca,'fontsize',fontsize)
                axis equal
                axis tight
            end

            % Shifting the spectrum such the center frequency is at the top left corner
            data_freq_complex_new_band(:,end) = [];
            data_freq_complex_new_band(end,:) = [];
            data_freq_complex_new_band = ifftshift(data_freq_complex_new_band);

            % Inverting back to the spatial domain
            data_band = ifft2(data_freq_complex_new_band,data_rows_new,data_columns_new);

            % Removing the padded region introduced when putting the number of lines and rows 2^n
            if data_rows-data_rows_new<0			% for the rows
                data_band(data_rows+1:end,:)=[];
            end
            if data_columns-data_columns_new<0		% for the columns
                data_band(:,data_columns+1:end)=[];
            end

            % Keeping only the real part
            data_band = real(data_band);

            % removing any mirrored image from the data
            if mirror_flag == 1
                data_band_or = data_band;
                data_band([1:n_mirror_y],:)=[];
                data_band([end-n_mirror_y+1:end],:)=[];
                data_band(:,[1:n_mirror_x])=[];
                data_band(:,[end-n_mirror_x+1:end])=[];
            end


            % plotting the final results when requested
            if plot_fig_flag==1
                h1= figure('name','output image');
                imagesc(X_lims_fig,Y_lims_fig,data_band)  
                % limit the colorbar to 95% bounds
                data_sorted = sort(reshape(data_band,[],1));
                ix = round([0.025 0.975].*length(data_sorted));
                colorlimits = [data_sorted(ix)];    
                caxis(colorlimits)
                colorbar
                set(gca,'fontsize',fontsize)
                xlabel('Distance [m]','fontsize',fontsize)
                ylabel('Distance [m]','fontsize',fontsize)
                title(['Band filter: ' , num2str(spatial_bands(kk,1)) ,' -- ' ,num2str(spatial_bands(kk,2)),' m'],'fontsize',fontsize)
                axis equal
                axis tight
                axis xy

                % saving of the figure when requested
                if save_fig==1
                    fig_save_name = [save_folder_str filesep 'output_' num2str(spatial_bands(kk,1)) '_'  num2str(spatial_bands(kk,2)) 'm.eps'];
                    set(h1,'PaperPositionMode','auto')
                    print(h1,'-depsc','-r150',fig_save_name)
                    clear fig_save_name
                    close(h1)
                end
                clear h1 



            end

            % Storing the output data 
            data_band_out(:,:,kk) = data_band;
            clear data_band
        end



        if k<=h_ifg_number && h_ifg_number~=1
            % thse are interferograms
            save_name = ['bandfilter_regular_hgt_ifg_' num2str(k) '.mat'];
        elseif k<=h_ifg_number && h_ifg_number==1
            % these are the heights
            save_name = 'bandfilter_regular_hgt.mat';
        else
            % thse are interferograms
             save_name = ['bandfilter_regular_ifg_' num2str(k-h_ifg_number) '.mat'];
        end


        % saving the bandfiltered data for eacht dataset seperately
        save([save_path, filesep, save_name],'data_band_out','x_res','y_res','spatial_bands','n_degree_butterworth','norm_filter_flag','dimension_filter')
        clear data_band_out 
    end
end





