function [P,h,T,DWPT,RH] = load_sounding_download(filename)
% [P,h,T,DWPT,RH] = load_sounding_profile_download(filename)
% Function to read and save the sounding data as donloaded by the sounding_profile_download function.
% Main program variables are saved as P [hPa], T [degrees], h [m] and  RH [%].
% The output file is saved under the same filename excluding the station number, i.e.
% it just contains the data YYYY_MM_DD_HH or acquisition.
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
% Bekaert David generated using the auto load function in matlab
% University of Leeds
%
% Modifications:
% DB 	03/2013 	Make comaptible with other sounding scripts
% DB	03/2013		Save the data automatically at the same path
% DB    10/2013   	Change filename 

startRow = 5;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray));
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
P = cell2mat(raw(:, 1));
h = cell2mat(raw(:, 2));
T = cell2mat(raw(:, 3));
DWPT = cell2mat(raw(:, 4));
RH = cell2mat(raw(:, 5));
MIXR = cell2mat(raw(:, 6));
DRCT = cell2mat(raw(:, 7));
SKNT = cell2mat(raw(:, 8));
THTA = cell2mat(raw(:, 9));
THTE = cell2mat(raw(:, 10));
THTV = cell2mat(raw(:, 11));


% saving the data in a mat variable. Omit the station number
[pathname,filename_temp,ext] = fileparts(filename);
pathname = [pathname,filesep];
clear ext
ix = find('_' == filename_temp);
if isempty(ix)==1
	fprintf('This is not the stationnumber_YYYYMMDD_HH.txt filename convention')
end
save_filename = [pathname filename_temp(ix(1)+1:end) '.mat'];
save(save_filename,'P','T','h','RH','MIXR','DRCT','SKNT','THTA','THTE','THTV')

%% Clear temporary variables
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
