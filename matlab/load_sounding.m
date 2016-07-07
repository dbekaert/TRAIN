function [] = load_sounding(sounding_data_path)
% Script which loads the downloaded sounding data into .mat files.
% Specify as input the full path to your sounding data e.g:
% >> sounding_data_path = '/nfs/a1/insar/italy/sounding_data/eu/16245';
% >> load_sounding(sounding_data_path)
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
% by David Bekaert -- University of Leeds
% modifications
% 21/10/2013 	DB 	Fixed change in file name of fuction, and suppress command line outputs
% 19/03/2014    DB  get the soundign datapath from the parm file

if nargin<1
    sounding_data_path = getparm_aps('sounding_dir');
end

% making a list with all the files that need to be saved into the matlab structure
external_command1 = ['echo files > ' sounding_data_path filesep 'file.list'];
external_command2 = ['ls ' sounding_data_path filesep '[0-9]*_????????_??.txt >> ' sounding_data_path filesep 'file.list'];

[a,b] = system(external_command1);
[a,b] = system(external_command2);
clear a b

% process all files specified by path and filename
files = char(textread([sounding_data_path filesep 'file.list'],'%s','headerlines',1));
n_files = size(files,1); 

% load each file in matlab and save it in the right structure
for k=1:n_files
	if floor(k/10)*10==k
		fprintf([num2str(k) ' completed out of ' num2str(n_files) '\n'])
	end
	load_sounding_download(files(k,:));
end
