function [matrix_data] = Dvec2Dmat(vector_data,WIDTH,LENGTH,IND,nodata)
% function which converts a vector to a matrix:
% option 1: WIDTH only where it is assumed that the original grid was 
%                  made using vector = reshape(matrix,[],1);
% option 2: WIDTH, LENGTH, IND where the IND has the same length as the vector_data
%                  and where IND relates to the positon in a matrix. Note
%                  that this method allows for a mask file to be included.
%                  You can use the matlab function IND = sub2ind([WIDTH LENGTH], I, J)
%                  to generate the IND variable
% nodata argument is optional and us used to fill the nodata. By default
% this is a NaN.
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
% modifications:
% 8/10/2016     DB  Add nodata option


%  checking of the input arguments
method =1;          % this is the width method
if nargin <2
    error(['Expecting two inputs: vector_data and WIDTH'])
end
if nargin>2 && nargin<4
    fprintf('WARNING: you need to specify 4 inputs when planning on using IND option, will fall back on width only \n')
end
if length(WIDTH)~=1
    error('WIDTH is given incorrect')
end
if nargin >=4
   if isempty(LENGTH) || isempty(IND) 
       error('You need to specify LENGTH and IND') 
   end
   if length(LENGTH)~=1
       error('LENGTH is given incorrect')
   end
   if length(IND)~=length(vector_data)
       error('Length of vector_data and IND needs to be consistent')
   end
   method = 2;      % change to the width, length and ind method
end
if nargin<5
    nodata=NaN;
end


%% storing the data again
if method ==1
    matrix_data = reshape(vector_data,WIDTH,[]);
elseif method ==2
    if isnan(nodata)
        matrix_data = NaN([WIDTH LENGTH]);
    elseif nodata==0
        matrix_data = zeros([WIDTH LENGTH]);
    else
        matrix_data = ones([WIDTH LENGTH]).*nodata;
    end
    matrix_data(IND)= vector_data;
end

