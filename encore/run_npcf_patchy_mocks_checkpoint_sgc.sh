#!/bin/bash

##################### DOCUMENTATION #####################
### Shell script for running the C++ encore NPCF-estimator function on a data and data-random catalog, then combining the outputs, including edge-correction (Oliver Philcox, 2021).
#
# This can be run either from the terminal or as a SLURM script (using the below parameters).
# The code should be compiled (with the relevant options, i.e. N-bins, ell-max and 4PCF/5PCF/6PCF) before this script is run. The isotropic 2PCF and 3PCF will always be computed.
# The script should be run from the code directory
# This is adapted from a similar script by Daniel Eisenstein.
# In the input directory, we expect compressed .gz files labelled {root}.data.gz, {ranroot}.ran.{IJ}.gz where {root} is a user-set name, and {IJ} indexes the random catalogs, from 0 - 31.
# We expect the summed weights to be the same for the data and each random catalog, but the random weights should be negative
# This script will compute the D^N counts, the (D-R)^N counts for 32 random subsets, and the R^N counts for one subset (should be sufficient).
# If the connected flag is set (and -DDISCONNECTED added to the makefile) we compute also the Gaussian contribution to the 4PCF.
# The output will be a set of .zeta_{N}pcf.txt files in the specified directory as well as a .tgz compressed directory of other intermediary outputs
# It is important to check the errlog file in the output directory to make sure all went well!
# Note that performing edge-correction is slow for the 5PCF and 6PCF since 9j symbols must be computed. Furthermore, the output multipoles are only accurate up to (ORDER-1), i.e. to compute an accurate edge-corrected spectrum with ell=5, we must compute (D-R) and R counts up to ell=6.
#
# NB: If needed, we could access a task ID by SLURM_ARRAY_TASK_ID, if we're running with SLURM
##########################################################
#SBATCH --array=1-100
#SBATCH -n 80 # cpus
#SBATCH -t 0-01:00:00 # time
#SBATCH --mem-per-cpu=1GB
#SBATCH -o /gpfs/akrolewski/parity_odd_4pcf/slurm_logs/%A.out         # File to which STDOUT will be written (make sure the directory exists!)
#SBATCH -e /gpfs/akrolewski/parity_odd_4pcf/slurm_logs/%A.err         # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification
#SBATCH --mail-user=alex@krolewski.com # Email to which notifications will be sent
#SBATCH -p preq

##################### INPUT PARAMETERS ###################

# Main inputs
useAVX=1 # whether we have AVX support
periodic=0 # whether to run with periodic boundary conditions (should also be set in Makefile)
connected=0 # if true, remove the disconnected (Gaussian) 4PCF contributions (need to set -DDISCONNECTED in the Makefile for this)
rmin=20 # minimum radius in Mpc/h
rmax=160 # maximum radius in Mpc/h

# Other inputs
scale=1 # rescaling for co-ordinates
ngrid=50 # grid-size for accelerating pair count
boxsize=1000 # only used if periodic=1

# File directories
ID=$(printf %04d $SLURM_ARRAY_TASK_ID)
echo "ID"
echo $ID
echo "SLURM ARRAY TASK ID"
echo $SLURM_ARRAY_TASK_ID
root="Patchy-Mocks-DR12SGC-COMPSAM_V6C_${ID}_cart" # root for data filenames
ranroot=Patchy-Mocks-DR12SGC-COMPSAM_V6C_x50_0001_cart # root for random filenames
ranroot0="Patchy-Mocks-DR12SGC-COMPSAM_V6C_x50_${ID}_cart"
in=/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/ # input directory (see above for required contents)
out="/gpfs/akrolewski/parity_odd_4pcf/out/10_bins/patchy_mocks/S/${SLURM_ARRAY_TASK_ID}/" # output file directory
tmp="/gpfs/akrolewski/scratch/patchy_mocks/S/${ID}" # temporary directory for intermediate file storage for this run (ideally somewhere with fast I/O)

# Load some python environment with numpy and sympy installed
module load miniconda3
source activate py4encore

##########################################################

# Define output coupling file (to avoid recomputation of disconnected pieces if multiple aperiodic simulations are run)
RRR_coupling=$out/$ranroot.RRR_coupling.npy

# Set number of threads (no SLURM)
#set OMP_NUM_THREADS = 4

# Set number of threads (with SLURM)
export OMP_NUM_THREADS=$SLURM_NPROCS

# Define command to run the C++ code
if (($useAVX == 1)); then
  code=./encoreAVX
else
  code=./encore
fi
echo $useAVX
echo $code

if (($periodic == 1)); then
  command="$code -rmax $rmax -rmin $rmin -ngrid $ngrid  -scale $scale -boxsize $boxsize"
else
  command="$code -rmax $rmax -rmin $rmin -ngrid $ngrid -scale $scale"
fi

mkdir $out

# Create a temporary directory for saving
#/bin/rm -rf $tmp       # Delete, just in case we have crud from a previous run.
export lineflag=0
if [ -f $out/progress.txt ]; then
    echo "progress.txt exists"
    while read -r line;
    do export lineflag=$line;
       done < $out/progress.txt
fi
echo "lineflag after read"
echo $lineflag

mkdir /gpfs/akrolewski/scratch/
mkdir $tmp

if (($lineflag == 0)); then
    
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi
if (($lineflag == 1)); then
    # Copy this script in for posterity
    cp $0 $tmp
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

# Create output directory
#if (!(-e $out)) then
#mkdir $out
#fi

errfile=errlog
errlog=$out/$errfile
tmpout=$tmp


# Create an output filefor errors
if (($lineflag == 2)); then
    #/bin/rm -rf $tmp
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

if (($lineflag == 3)); then
    date > $errlog
    echo Executing $0 >> $errlog
    echo $command >> $errlog
    echo $OMP_NUM_THREADS >> $errlog
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

# Filename for saved multipoles (a big file)
multfile=$tmp/$root.mult

# Extra the data into our temporary ramdisk
#gunzip -c $in/$root.data.gz > $tmp/$root.data

# Find number of galaxies (needed later for R^N periodic counts)
Ngal=`cat $in/$root.data | wc -l`
Ngal=`expr $Ngal + 1`
echo $lineflag
echo "Lineflag at 143"

#### Compute D^N NPCF counts
# Note that we save the a_lm multipoles from the data here
if (($lineflag == 4)); then
    echo Starting Computation
    echo "Starting D^N" >> $errlog
    date >> $errlog
    echo "($command -in $in/$root.data -save $multfile -outstr $root.data > $tmpout/$root.d.out) &>> $errlog"
    ($command -in $in/$root.data -save $multfile -outstr $root.data > $tmpout/$root.d.out) &>> $errlog
    echo "Done with D^N"
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

# Remove the output - we don't use it
#rm output/$root.data_?pc*.txt



### Compute R^N NPCF counts
# We just use one R catalog for this and invert it such that the galaxies are positively weighted
#gunzip -c $in/$ranroot.ran.00.gz > $tmp/$root.ran.00
if (($lineflag == 5)); then
    #cp $in/$ranroot.ran.00 $tmp/$root.ran.00
    cp ${in}/$ranroot0.ran.00 $tmp/$root.ran.00
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

if (($lineflag == 6)); then
    echo "Starting R^N" >> $errlog
    date >> $errlog
    echo "($command -in $tmp/$root.ran.00 -outstr $root.r -invert > $tmpout/$root.r.out) &>> $errlog"
    ($command -in $tmp/$root.ran.00 -outstr $root.r -invert > $tmpout/$root.r.out) &>> $errlog
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi
# Copy the output into the temporary directory
if (($lineflag == 7)); then
    mv output/$root.r_?pc*.txt $tmpout/
    echo "Done with R^N"
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi



# Now make D-R for each of 32 random catalogs, with loading
for n in  00 01 02 03 04 05 06 07 08 09 \
	    10 11 12 13 14 15 16 17 18 19 \
	    20 21 22 23 24 25 26 27 28 29 \
	    30 31 
do
    # First copy the randoms and add the data
    #/bin/cp -f $tmp/$root.data $tmp/$root.ran.$n
    #gunzip -c $in/$ranroot.ran.$n.gz >> $tmp/$root.ran.$n
    echo "Line 207" >> $errlog
    echo $lineflag >> $errlog
    if (($lineflag == $((10#$n*3 + 8)))); then
        cat "$in/$root.data" "$in/$ranroot.ran.${n}" > "$tmp/$root.ran.${n}"
        echo $lineflag >> $out/progress.txt
        export lineflag=$(($lineflag + 1))
    fi

    echo "Line 215" >> $errlog
    echo $lineflag >> $errlog

    ### Compute the (D-R)^N counts
    # This uses the loaded data multipoles from the D^N step
    # Note that we balance the weights here to ensure that Sum(D-R) = 0 exactly
    if (($lineflag == $((10#$n*3 + 9)))); then
        echo "Starting D-R ${n}" >> $errlog
        date >> $errlog
        echo "($command -in $tmp/$root.ran.${n} -load $multfile -outstr $root.n${n} -balance > $tmpout/$root.n${n}.out)" &>> $errlog
        ($command -in "$tmp/$root.ran.${n}" -load $multfile -outstr "$root.n${n}" -balance > "$tmpout/$root.n${n}.out") &>> $errlog
        echo $lineflag >> $out/progress.txt
        export lineflag=$(($lineflag + 1))
    fi
    echo "Line 229" >> $errlog
    echo $lineflag >> $errlog

    # Copy the output into the temporary directory
    if (($lineflag == $((10#$n*3 + 10)))); then
        mv output/$root.n${n}_?pc*.txt $tmpout/
        echo "Done with D-R ${n}" >> $errlog
        echo $lineflag >> $out/progress.txt
        export lineflag=$(($lineflag + 1))
    fi

    # Remove the random catalog
    #/bin/rm -f $tmp/$root.ran.$n
    

done # foreach D-R loop

### Now need to combine the files to get the full NPCF estimate
# We do this in Python, and perform edge-correction unless the periodic flag is not set
if (($lineflag == 104)); then
    if (($periodic == 1)); then
      echo 'Combining files together without performing edge-corrections (using analytic R^N counts)'
      python python/combine_files_periodic.py $tmpout/$root $Ngal $boxsize $rmin $rmax &>> $errlog
    else
      echo 'Combining files together and performing edge-corrections'
      python python/combine_files.py $tmpout/$root &>> $errlog $OMP_NUM_THREADS
    fi
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

### If the connected flag is set, also combine files to estimate the disconnected 4PCF
# We perform edge corrections unless the periodic flag is not set.
# If the file RRR_coupling exists we load the edge-correction matrix from file, else it is recomputed
if (($lineflag == 105)); then
    if (($connected == 1)); then
      if (($periodic == 1)); then
        echo 'Combining files together to compute the disconnected 4PCF without performing edge corrections'
        python python/combine_disconnected_files_periodic.py $tmpout/$root 4 $Ngal $boxsize $rmin $rmax &>> $errlog
      else
        echo 'Combining files together to compute the disconnected 4PCF including edge corrections'
        python python/combine_disconnected_files.py $tmpout/$root 4 $RRR_coupling &>> $errlog
      fi
    fi
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

# Do some cleanup
#rm $multfile

# Now move the output files into the output directory.
# Compress all the auxilliary files and copy
echo 'Finished with computation.  Placing results into $out/'
echo 'Finished with computation.  Placing results into $out/ >> $errlog'
if (($lineflag == 106)); then
    date >> $errlog
    pushd $tmpout > /dev/null
    echo >> $errlog
    /bin/ls -l >> $errlog
    /bin/cp $errlog .
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi
if (($lineflag == 107)); then
    tar cfz $root.tgz $root.*.out $root.*pc*.txt $errfile $0
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi
if (($lineflag == 108)); then
    popd > /dev/null
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi
if (($lineflag == 109)); then
    /bin/mv $tmpout/$root.tgz $tmpout/$root.zeta_*pcf.txt $out/
    #/bin/rm -rf $tmp
    echo $lineflag >> $out/progress.txt
    export lineflag=$(($lineflag + 1))
fi

# Destroy ramdisk
#

