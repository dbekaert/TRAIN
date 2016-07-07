function [y_mean,y_std_mean,obs_kept] = w_test_mean(y,y_std,alpha)
% [y_mean,y,y_std,obs_kept] = w_test_mean(y,y_std,alpha)
% Function which performs an outlier detection for the mean of a set of
% observations y and corresponding standard deviations y_std. The rate of
% rejection is determined by alpha. Alpha represents the probability that a
% good observations was considered as an outlier, i.e. alpha=1 all is rejected. 
%
% inputs: 
% y:        Observations, this needs to be a column or row vector of length n.
% y_std:    standard deviations of the observations. Needs to be a column
%           or row vector with length n. By default equal weighting is
%           assumed.
% alpha:    The probability that a point while good was regarded as an
%           outlier. When 1 all will be rejected, the smaller a less
%           outlier rejection criteria is set.
% 
% outputs:
% y_mean:   The weighted mean after rejection of the outliers.
% obs_kept: Vector with the position of those observations kept
%
%     w_test_mean.m
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
% By David Bekaert - February 2013
% PhD student - University of Leeds


if nargin<1
    
end
if size(y,2)>1
    y = y';
    if size(y,2)>1
        error('myApp:argChk', ['number of observastions y and corresponding standard deviation are inconsistent...\n'])   
    end
end
if nargin<2 || isempty(y_std)
    y_std = ones(size(y));
end
if nargin<3 || isempty(alpha)
   alpha = 0.05;  
end


% computation of the rejection threshold k_alpha_thres. 
% this is based on a standard normal distribution with level of
% significance alpha. alpha represents the probability that a point 
% is regarded as outlier while in reality it was not.
k_alpha_thres = norminv(1-alpha,0,1);



% a variable to trace which observations were rejected
obs_kept = 1:length(y);

continue_reject = 1;
while continue_reject
    % Total number of observations to start with
    n_obs = length(y);

    % Variance-matrix of the null-hypothesis H0
    Qyy = diag(y_std).^2;

    % The design matrix of the model (the mean)
    A = ones([n_obs 1]);
    
    % inverse variance matrix
    Qyy_inv = inv(Qyy);
    
    % estimate the mean
    x_est = (A'*Qyy_inv*A)\A'*Qyy_inv*y;
        
    % comptuation of the observation residuals
    y_est = A*x_est;
    e_est = y - y_est;
    
    % propagation of the variances  
    Qx_estx_est= inv(A'*Qyy_inv*A);         % variance matrix of the estimated unknowns
    Qy_esty_est = A*Qx_estx_est*A';         % variance matrix estimated observations
    Qe_este_est = Qyy - Qy_esty_est;        % variance matrix residuals
    
    % computation of the weights
    w = NaN([n_obs 1]);
    for k=1:n_obs
        % setting the cy colum vector
        cy = zeros([n_obs 1]);
        cy(k)=1;
        
        % compuation of the corresponding weight
        w(k)=(cy'*Qyy_inv*e_est)./(sqrt(cy'*Qyy_inv*Qe_este_est*Qyy_inv*cy));
        clear cy
    end
    
    
    % take the absolute w-weight
    w = abs(w);
    
    % Find the biggest w value and test it for the hypothesis 
    ix = find(max(w)==w);
    if length(ix)>1
        ix=ix(1);
        % this is not ideal and needs revising
    end
    
    
    % Check if the hypothesis is satisfield 
    if w(ix)> sqrt(k_alpha_thres)         
       
        % the test rejected all observations. The variance model is likely
        % to be incorect
        if n_obs==1
            continue_reject=0;
        else
            % removing the identified observation as it appears and outlier
            y(ix) = [];         % update the observations
            y_std(ix,:)=[];     % update the design matrix
            obs_kept(ix)=[];    % remove the outlier from the list
            clear ix
        end
    else
        clear ix
        % no more outliers found, end the outlier rejection
        continue_reject=0;
    end

end

% Once all outlier have been rejected y_est and y_est_std give the left
% over dataset. The weighted mean is represented by x_est
y_mean = x_est;
% The standard deviation of the weighted mean is contained in Qx_estx_est
y_std_mean = sqrt(Qx_estx_est);












