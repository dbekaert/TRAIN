function [parm_success] = get_parm_xml(xmlfile,parm,parm_key)
% function which will read xml file and updates the parameter_field matching
% the requested parameter. Will return "FAIL" or "SUCCESS" based on
% whether it found the variable
%
% Bekaert David - Jet Propulsion Laboratory

if nargin<3
    error('Require three inputs')
end

% initialize the field
temp = 0;
count_changes = 0;

% loop over the xml file 
fid = fopen(xmlfile,'r');
fid_new = fopen([xmlfile '_temp'],'w');


while temp~=-1
    temp = fgetl(fid);
    if temp==-1
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
            % write the output in case its not a value line
            if isempty(ix_value);
                fprintf(fid_new,[temp '\n']);
            end
            
            temp = fgetl(fid);
            ix_value = strfind(temp,'<value>');
            
            % avoid inf loop
            if temp==-1
                break
            end
        end

        if temp~=-1
            % update with new value or string
            if isnumeric(parm_key)
                parm_key = strtrim(num2str(parm_key));
            end
            parm_key_str = [repmat(' ',1,ix_value-1) '<value>' parm_key '</value>'];
            fprintf(fid_new,[parm_key_str '\n']);
            count_changes = count_changes+1;
        end
    else
        fprintf(fid_new,[temp '\n']);
    end
end
fclose(fid);
fclose(fid_new);

% rename the file
delete(xmlfile)
movefile([xmlfile '_temp'],xmlfile)

% see if the search was succesful
if count_changes==0
    fprintf(['Could not find ' parm '... \n'])
else
    parm_success = 'SUCCESS';
end