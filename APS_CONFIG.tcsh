#     Copyright (C) 2015  Bekaert David - University of Leeds
#     Email: eedpsb@leeds.ac.uk or davidbekaert.com
# 
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#


# Give the correct path to the APS toolbox
setenv APS_toolbox       "/Users/dbekaert/Software/svn_aps/trunk/"
setenv PYTHONPATH	 "/Users/dbekaert/Software/svn_aps/trunk/python_modules/"
# full path to the get_modis.py file 
setenv get_modis_filepath	"/Users/dbekaert/Software/oscar/get_modis.py" 	
setenv get_ecmwf_filepath 	"/Users/dbekaert/Software/ecmwfapi/api.py"

# when matlab search for you shell.
# if not found a default bash shell is used. 
# uncomment the next line and specify your shell when needed. 
setenv MATLAB_SHELL  "/bin/csh"

#####################################
# shouldn't need to change below here
#####################################

if (! ${?MATLABPATH}) then
        setenv MATLABPATH "$APS_toolbox/matlab"
else if (":${MATLABPATH}:" !~ *":$APS_toolbox/matlab:"*) then
        setenv MATLABPATH "$APS_toolbox/matlab:${MATLABPATH}"
endif

setenv APS_toolbox_scripts "$APS_toolbox/scripts"
setenv APS_toolbox_bin  "$APS_toolbox/bin"
setenv PATH "${PATH}:${APS_toolbox_bin}:${PYTHONPATH}"


