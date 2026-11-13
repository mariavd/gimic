#!/bin/bash

SCRIPTS_DIR=$(pwd)

pathToGimic=$(which gimic)
if [ -z $pathToGimic ]
then
    echo "The path to the GIMIC binary is not defined. Please add it to your .bashrc as"
    echo export PATH=\"/path/to/gimic/\":\${PATH}
fi


# The current profile scripts (current/vertical/squares-profile-*.sh, gimic-run.sh,
# plot-current-profile.sh, crit_pts.sh) are thin wrappers around gimic-profile,
# which is installed together with gimic; nothing to prepare for them.

# Prepare all other scripts
SCRIPTS_IN=$(ls src/*in)

for file in $SCRIPTS_IN
do
    SCRIPT_OUT=$(echo $file | sed -e s:^src/:: -e s:.in:: )
    sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
    if [ $SCRIPT_OUT -ot $file ] 
    then 
        sed "$sedstring" $file > $SCRIPT_OUT
        echo "Created script $SCRIPT_OUT."
	chmod +x $SCRIPT_OUT
    fi
done

echo; echo "Path to GIMIC:"
echo "export GIMIC_HOME=\"$SCRIPTS_DIR\""
