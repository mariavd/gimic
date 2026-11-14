#!/bin/bash

SCRIPTS_DIR=$(pwd)

# NOTE: "make install" of GIMIC already configures and installs these scripts
# next to gimic. This script is only needed to use them from this directory.

pathToGimic=$(command -v gimic)
if [ -z "$pathToGimic" ]
then
    echo "The path to the GIMIC binary is not defined. Please add it to your .bashrc as"
    echo export PATH=\"/path/to/gimic/bin\":\${PATH}
fi

GIMIC_PREPARE=$(command -v gimic-prepare)
if [ -z "$GIMIC_PREPARE" ]
then
    echo "gimic-prepare is not on the PATH: the profile scripts will look for it at run time."
    GIMIC_PREPARE=gimic-prepare
fi

# Prepare all other scripts
SCRIPTS_IN=$(ls src/*in)

for file in $SCRIPTS_IN
do
    SCRIPT_OUT=$(echo $file | sed -e s:^src/:: -e s:.in:: )
    sedstring="s:@GIMIC_PREPARE@:$GIMIC_PREPARE:"
    if [ $SCRIPT_OUT -ot $file ] 
    then 
        sed "$sedstring" $file > $SCRIPT_OUT
        echo "Created script $SCRIPT_OUT."
	chmod +x $SCRIPT_OUT
    fi
done

echo; echo "Path to GIMIC:"
echo "export GIMIC_HOME=\"$SCRIPTS_DIR\""
