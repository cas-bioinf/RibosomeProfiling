#!/bin/sed 6,7!d;s/^# //;s/ *#$//

########################################################################################################################
# These common variables are expected to be set before running of the other scripts.                                   #
#                                                                                                                      #
# This script cannot be run directly ('./0-variables.sh') as that does not set variables outside the script.           #
# Run the script via '. ./0-variables.sh' or 'source ./0-variables.sh', or copy all commands manually.                 #
#                                                                                                                      #
# Paths should end with a directory separator ('/') as variables are used as "${VAR}filename" (instead of              #
# $VAR/filename) so if the files are in the working directory, it is not necessary to define $VAR.                     #
#                                                                                                                      #
# It is recommended to use absolute paths for the variables are valid even if the working directory is changed         #
# ('realpath' command resolves absolute paths, however the relative paths must still be valid during definition)       #
########################################################################################################################

# Location of all programs that are not installed globally
programs=$( realpath programs )/
echo '$programs   = '"'$programs'" 

# Location of reference sequences like rRNA, tRNA, genome etc.
references=$( realpath references )/
echo '$references = '"'$references'"

# Location of raw sequencing data
input=$( realpath input )/
echo '$input      = '"'$input'"

# Where to store output files
output=$( realpath output )/
echo '$output     = '"'$output'"

# Where to store logs
logs=$( realpath log )/
echo '$logs       = '"'$logs'"
