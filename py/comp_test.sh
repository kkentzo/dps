#!/bin/bash

# USAGE:
# provide check argument to just check whether all files 
# were created 

#SBATCH --job-name=COMP
#SBATCH --output=comp.out
#SBATCH --error=comp.err
#SBATCH --ntasks=100



# prefix arguments to DPS executable
LD_LIBRARY_PATH=/home/kkentzo/local/lib
# the path to dps
DPS_EXE=../src/dps
# the path to the input HDF files (containing
# the constructed populations)
INPUT_DIR=./populations/c
# the path to the results files
RESULTS_DIR=./results
# the value of PCONJ
PCONJ=0.01
# common arguments to dps
ARGS="-c --steps 100 --print_every 0"
# how many cores to use
CORES=4
# how may runs per value of pconj
RUNS=1

# the number of KAPPA and ALPHA values
N_KAPPA=10
N_ALPHA=10

echo "======================================="
echo "SAVING TO : $RESULTS_DIR"
echo "LOADING FROM : $INPUT_DIR"

echo "KAPPA_VALUES=$N_KAPPA"
echo "ALPHA_VALUES=$N_ALPHA"
echo "PCONJ=$PCONJ"
echo "RUNS=$RUNS"
echo "CORES=$CORES"
echo "TOTAL NUMBER OF JOBS=`expr $N_KAPPA \* $N_ALPHA \* $RUNS`"
echo "======================================="

# stamp the date
date

# count the jobs launched so far (for properly inserting waits)
JOB_COUNTER=0

# create results directory
mkdir -p $RESULTS_DIR

# ==============================================================
# VERIFY THAT ALL THE INPUT FILES EXIST
for idx_kappa in $(seq 0 $((${N_KAPPA} - 1)))
do
    for idx_alpha in $(seq 0 $((${N_ALPHA} - 1)))
    do
	# determine input file name
	IN_FNAME=$INPUT_DIR/population.$idx_kappa.$idx_alpha.h5
	if [ ! -f $IN_FNAME ]
	then
	    echo "Input file does not exist! : $IN_FNAME"
	    echo "Aborting."
	    exit
	fi
    done
done

# ==============================================================
# SCHEDULE SIMULATIONS
for idx_kappa in $(seq 0 $((${N_KAPPA} - 1)))
do
    for idx_alpha in $(seq 0 $((${N_ALPHA} - 1)))
    do
	# determine input file name
	IN_FNAME=$INPUT_DIR/population.$idx_kappa.$idx_alpha.h5

	# conduct all runs
	for idx_run in $(seq 0 $(($RUNS-1)));
	do
	    # determine output file name
	    OUT_FNAME=$RESULTS_DIR/results.$idx_kappa.$idx_alpha.$idx_run.h5

	    # SCHEDULE SIMULATIONS
	    if [ `expr $JOB_COUNTER % $CORES` -eq 0 ]
	    then
		# insert a wait 
		wait
	    fi
	    
	    if [ ! -f $OUT_FNAME ]
	    then
	        # launch process
		LD_LIBRARY_PATH=$LD_LIBRARY_PATH $DPS_EXE $ARGS --pconj $PCONJ \
		    --load_from $IN_FNAME $OUT_FNAME &
	
	        # increment job counter
		JOB_COUNTER=`expr $JOB_COUNTER + 1`
	    else
		# file exists - do **not** overwrite
		echo "SKIPPING FILE $OUT_FNAME"
	    fi
	done
    done
done

# insert a final wait (for the last batch of jobs)
wait

# stamp the date
date
