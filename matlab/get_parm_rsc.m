function [parm_key,parm_key_comp] = get_parm_rsc(rscfile,parm)
% function which will read rsc file and find the parameter_field matching
% the requested parameter. Will return string and numbers depending on the
% parameter_field
%
% Bekaert David - Jet Propulsion Laboratory
% modifications:
% DB    8/10/2016   Make sure not to end into an inf loop for the second while statement
% DB    3/11/2017   Trying to cpature incase there are mutiple variables with same parm matching

if nargin<2
    error('Require two inputs')
end

% initialize the field
counter=1;

% loop over the xml file 
fid = fopen(rscfile,'r');
temp=1;
parm_key = 'FAIL';


while temp~=-1
    temp = fgetl(fid);
    if temp==-1
        break
    end
    
    % search for the variable
    ix =  strfind(temp,[parm]);
    if isempty(ix)
        ix = strfind(temp,[lower(parm)]);
    end
    if isempty(ix)
        ix = strfind(temp,[upper(parm)]);
    end
    if ~isempty(ix)
        % see if the value is at this string or not
        ix_value = length(parm);
        % remove the start and trim the remaining
        temp(1:ix_value)=[];
        parm_key = strtrim(temp);
        break
    end
end
fclose(fid);



% see if the search was succesful
if strcmpi(parm_key,'FAIL')
    fprintf(['Could not find ' parm '... in .rsc file\n'])
else
    if ~isempty(str2num(parm_key))
        parm_key= str2num(parm_key);
    end
end