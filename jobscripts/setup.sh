#!/bin/bash

SCRIPTS_DIR=$(pwd)

checkMaxProj=$(command -v maximise_projection)
if [ ! -e $checkMaxProj ]
then
    echo "The program maximise_projection is not found."
    echo "Compile it and export its path to use the automatic orientation of the magnetic field."
#    exit
fi

pathToGimic=$(which gimic)
if [ -z $pathToGimic ]
then
    echo "The path to the GIMIC binary is not defined. Please add it to your .bashrc as"
    echo export PATH=\"/path/to/gimic/\":\${PATH}
fi


# Prepare the profile scripts: a header (which parameters to ask for) and a
# submit part (run locally with OpenMP threads, or submit to SLURM).
for mode in current vertical squares
do
    for where in local cluster
    do
        file1="src/$mode-profile-header"
        file2="src/current-profile-$where-submit"
        SCRIPT_OUT="src/$mode-profile-$where.sh.in"
        if [ $SCRIPT_OUT -ot $file1 ] || [ $SCRIPT_OUT -ot $file2 ] || [ $SCRIPT_OUT -ot src/current-profile-header ] || [ $SCRIPT_OUT -ot src/functions-def ]
        then
            cat $file1 > $SCRIPT_OUT
            if [ $mode != current ]
            then
                cat src/current-profile-header >> $SCRIPT_OUT   # the vertical and squares profiles are the same script with another slicing
            fi
            cat $file2 >> $SCRIPT_OUT
        fi
    done
done

# Prepare the batch job scripts:

file="src/jobscript.IN"
SCRIPT_OUT=jobscript
sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
if [ $SCRIPT_OUT -ot $file ] || [ $SCRIPT_OUT -ot src/jobscript-header ]
then 
    echo; echo "REMEMBER TO CHANGE THE BATCH SCRIPT jobscript-header TO SUIT YOUR CLUSTER BEFORE SETUP"
    cat src/jobscript-header > $SCRIPT_OUT
    sed "$sedstring" $file >> $SCRIPT_OUT
    echo "Created script $SCRIPT_OUT."
    chmod +x $SCRIPT_OUT
fi

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
