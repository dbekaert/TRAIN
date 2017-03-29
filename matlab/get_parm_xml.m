function [parm_key] = get_parm_xml(xmlfile,parm)
% function which will read Rxml file and find the parameter_field matching
% the requested parameter. Will return string and numbers depending on the
% parameter_field
%
% Bekaert David - Jet Propulsion Laboratory
% modifications:
% DB    8/10/2016   Make sure not to end into an inf loop for the second while statement

if nargin<2
    error('Require two inputs')
end

% initialize the field
parm_key = [];

% loop over the xml file 
fid = fopen(xmlfile,'r');
while isempty(parm_key) 
    temp = fgetl(fid);
    if temp==-1
        parm_key = 'FAIL';
        break
    end
    ix =  strfind(temp,['"' parm '"']);
    if isempty(ix)
        ix = strfind(temp,['"' lower(parm) '"']);
    end
    if ~isempty(ix)
        % see if the value is at this string or not
        ix_value = strfind(temp,'<value>');
        while isempty(ix_value);
            temp = fgetl(fid);
            ix_value = strfind(temp,'<value>');
            if temp==-1
                parm_key = 'FAIL';
                break
            end
        end
        
        if temp~=-1
            % remove the start of <value>
            temp(1:ix_value+length('<value>')-1)=[];
            % remove the end of value
            ix_value = strfind(temp,'</value>');
            temp(ix_value:end)=[];
            parm_key = strtrim(temp);
        end
    end
end
fclose(fid);


% see if the search was succesful
if strcmpi(parm_key,'FAIL')
    fprintf(['Could not find ' parm '... \n'])
else
    % try to see if this is a number of not
    if ~isempty(str2num(parm_key))
        parm_key = str2num(parm_key);
    end
end