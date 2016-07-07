function []=setparm_aps(parmname,value,newflag)
%SETPARM_APS sets parameters in a saved workspace
% SETPARM_APS(PARMNAME,VALUE,FLAG) 
% Only enough characters of PARMNAME to make it unique need be typed.
% If VALUE is set to nan, the parameter is reset to the default value.   
% FLAG is optional, valid values are:
%       1 = add a new parameter (PARMNAME must be typed in full)
%      -1 = delete parameter from workspace(VALUE ignored)
%
% Based on script by Andy Hooper (StaMPS)
% Modified for the aps toolbox by David Bekaert - University of Leeds - 2013
% This script allows for non-StaMPS structured processed data.
%
% modifications:
% 04/2013   DB:     Convert column values to a row when logging the changes
% 04/2013   DB      Allowing the grid resolution to be specified by single
%                   value.
% 02/2015   DB      Fix in case you give only H:MM for UTC_sat
% 11/2015   DB      Include more freedom in specifying UTC_sat
% 02/2016   DB 	    Give better output for the band selection

parmfile='parms_aps.mat';
localparmfile='localparms_aps.mat';

parent_flag=0; % 0 when parms.mat in current directory


if exist(['.' filesep parmfile],'file')==0 && exist(['..' filesep parmfile],'file')==0
    parms_default_aps
end
if exist(['.' filesep parmfile],'file')
    parms=load(parmfile);
elseif exist(['..' filesep parmfile],'file')
    parmfile=['..' filesep parmfile];
    parms = load(parmfile);
    parent_flag=1; % 1 when parms.mat in parent directory
end

if exist(localparmfile,'file')
    localparms=load(localparmfile);
else
    localparms=struct('Created',date);
end



% a flag which changes the band of the powerlaw and updates the data file
% as well
if nargin>2 & newflag==1
    if isempty(value)
        logit([parmname,' = []'],'aps.log',parent_flag);
    elseif isnumeric(value)
        logit([parmname,' = ',num2str(value)],'aps.log',parent_flag);
    else
        logit([parmname,' = ',value],'aps.log',parent_flag);
    end
    parms=setfield(parms,parmname,value);
    save(parmfile,'-struct','parms')
      
else
    if nargin>1
        if ~isnumeric(parmname)
            parmnum=strmatch(parmname,fieldnames(parms)); 
            if length(parmnum)>1
                error(['Parameter ',parmname,'* is not unique'])
            elseif isempty(parmnum)
                error(['Parameter ',parmname,'* does not exist'])
            end
        else
            parmnum=parmname;
        end
        parmnames=fieldnames(parms);
        if size(parmnames,1)<parmnum
            error(['There are only ',num2str(size(parmnames,1)),' fields'])
        end
        parmname=parmnames{parmnum};
        
        if isnan(value)
            parms=rmfield(parms,parmname);
            save(parmfile,'-struct','parms')

            if isfield(localparms,parmname)
                localparms=rmfield(localparms,parmname);
                save(localparmfile,'-struct','localparms')
            end

            parms_default_aps
            value=getparm_aps(parmname);
               
            % updating the data variable
            if strcmpi(parmname,'powerlaw_kept')
                
                if ~strcmpi(getparm_aps('powerlaw_all_bands'),'y')
 		   fprintf('There are no bands stored, keep the original \n')
                else
                    aps_powerlaw_update_band(value)
                end
            end
            disp([parmname,' reset to default value'])
            

            
        else
            if nargin<=2 || newflag>=0
              if isempty(value)
                    logit([parmname,' = []'],'aps.log',parent_flag); 
              elseif isnumeric(value)  
                   if strcmp(parmname,'powerlaw_xy_res') && size(value,2)==1
                        value = [value value];
                   end
                   % updating the data variable
                   if strcmpi(parmname,'powerlaw_kept')
                        if ~strcmpi(getparm_aps('powerlaw_all_bands'),'y')
 				fprintf('There are no bands stored, keep the original \n')                            
                        else
                            aps_powerlaw_update_band(value)
                        end
                   end
                   
                   if size(value,1)>1                   
                       value_str = '[';
                       for ll=1:size(value,1)
                           if ll==size(value,1)
                               value_str = [value_str num2str(value(ll,:)) ];
                           else
                                value_str = [value_str num2str(value(ll,:)) '; '];
                           end
                       end
                       value_str = [value_str ']'];
                   else
                       value_str = num2str(value);
                   end
                   logit([parmname,' = ',value_str],'aps.log',parent_flag);
                   clear value_str                           
              else
                  if size(value,1)>1
                      value_str = '[';
                      for ll=1:size(value,1)
                           if ll==size(value,1)
                               value_str = [value_str num2str(value(ll,:)) ];
                           else
                                value_str = [value_str num2str(value(ll,:)) '; '];
                           end
                       end
                       value_str = [value_str ']'];
                  else
                       if strcmpi(parmname,'UTC_sat') & (size(value,2)~=5 & size(value,2)~=11)
                          fprintf('UTC_sat is not given as ''HH:MM'', fixing this...\n')
                          % the UTC time does not have HH:MM 
                          ix = find(value==':');
                          value_temp = '00:00';
                          value_temp(3-ix+1:3-ix+size(value,2))=value;
                          value = value_temp;
                          clear value_temp
                       end
                       value_str = value;
                   end
                   logit([parmname,' = ',value_str],'aps.log',parent_flag);
                   clear value_str
              end
            end
            
            if nargin>2 
                if newflag==2
                    localparms=setfield(localparms,parmname,value);
                    save(localparmfile,'-struct','localparms')
                    logit('Added to LOCAL parameter file','aps.log');
                elseif newflag==-1
                    parms=rmfield(parms,parmname);
                    save(parmfile,'-struct','parms')
                    logit([parmname,' removed from parameter file'],'aps.log',parent_flag);
                elseif newflag==-2
                    localparms=rmfield(localparms,parmname);
                    if size(fieldnames(localparms),1)>1
                        save(localparmfile,'-struct','localparms')
                    else
                        delete(localparmfile)
                    end
                    logit([parmname,' removed from LOCAL parameter file'],'aps.log');

                else
                    error('Invalid value for NEWFLAG')
                end
            else
                if ~isfield(localparms,parmname)
                    parms=setfield(parms,parmname,value);
                    save(parmfile,'-struct','parms')
                else
                    localparms=setfield(localparms,parmname,value);
                    save(localparmfile,'-struct','localparms')
                    disp('Warning: Only LOCAL parameter file updated')
                end
               
            end
        end
    elseif nargin >0
        error('Format is: SETPARM(PARMNAME,VALUE,[NEWFLAG])')
    else
        disp(orderfields(parms))
        if size(fieldnames(localparms),1)>1
            localparms
        end
    end
end

