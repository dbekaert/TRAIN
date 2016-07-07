function [error_flag] = aps_systemcall(command_str);
% Make a systemcall an check for errors too
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
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.%
%
% by Bekaert David
% modifications
% DB    02/2016     Add exeption error_flag to be passed up.

% call the system
[system_call_argument1,system_call_argument2]=system(command_str);


% check for the error message
[error_flag] = aps_check_systemcall_error(system_call_argument1,system_call_argument2);
