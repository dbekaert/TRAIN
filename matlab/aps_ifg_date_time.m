function [ifg_dates_out,UTC_sat] = aps_ifg_date_time(ifg_dates,UTC_sat)

% change definition of UTC time to allow for one for each date:
% user could do UTC1&UTC2 for each ifg, or given only 1 UTC
n_ifg = size(ifg_dates,1);
% redefine UTC if needed
if size(UTC_sat,1)==1 & size(UTC_sat,2)==5
    UTC_sat = repmat([UTC_sat '&' UTC_sat],n_ifg,1);
end

% Check if UTC was defined correctly.
[ix_row ix_column] = find(UTC_sat=='&');
if sum(abs(diff(ix_column)))~=0
    fprintf('Use convention "UTC1&UTC2" for each IFG on one line, or single line "UTC" in case same UTC time for all IFG');
    error('Make sure that the UTC_sat time has the & seperation') 
end

% consistency check
n_ifg_utc = size(UTC_sat,1);
if n_ifg_utc>1 & n_ifg_utc~=n_ifg
    fprintf('Use convention "UTC1&UTC2" for each ifg one line or "UTC" in case same utc time for all');
    error('Your specified more than 1 UTC time, but it did not match your number of IFGs\n')    
end

% add UTC time to the IFG dates
ifg_dates_out = [datenum([datestr(ifg_dates(:,1),'yyyymmdd') UTC_sat(:,[1 2 4 5])],'yyyymmddHHMM')   datenum([datestr(ifg_dates(:,2),'yyyymmdd') UTC_sat(:,[7 8 10 11])],'yyyymmddHHMM')];







