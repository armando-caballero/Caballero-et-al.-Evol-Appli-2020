#!/bin/bash
#$ -cwd

rm script_SLIM_ID_100.sh.*
rm data*
rm qt.phe
rm list*
rm outfileSLIM
rm timefile
rm checkfile
rm slimout 2> /dev/null

#Check number of arguments
if [ $# -ne 4 ]  
then
	echo "Usage: $0 <INPUT> <NIND> <REPS> <d>" 
	exit 1
fi

### Set arguments

INPUT=$1
NIND=$2
REPS=$3
n=$4

WDIR=$PWD 
mkdir -p /state/partition1/slimDruet$n/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################

cp shell_F_values /state/partition1/slimDruet$n/$SLURM_JOBID/
cp plink1.9 /state/partition1/slimDruet$n/$SLURM_JOBID/
cp gcta64 /state/partition1/slimDruet$n/$SLURM_JOBID/
cp slim-bivariate /state/partition1/slimDruet$n/$SLURM_JOBID/
cp seedfile /state/partition1/slimDruet$n/$SLURM_JOBID/
cp $INPUT /state/partition1/slimDruet$n/$SLURM_JOBID/
cp SNP_BP /state/partition1/slimDruet$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/slimDruet$n/$SLURM_JOBID

################## LOOP FOR REPLICATES ####################

module load gsl/2.1

for((i=1;i<=$REPS;i++))
do

START=$(date +%s)
./slim-bivariate $INPUT
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Slim took 		$DIFF seconds" > timefile

###################### SNP_BP #########################

START=$(date +%s)
#0:fitness, 1:QT; N; SDe
./SNP_BP<<@
0
-99
1	trait (0:fitness; 1: QT)
$NIND	N
0.18	scale_a
0.0	Ve
0.05 MAF for HML
@
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_BP took 	$DIFF seconds" >> timefile

cat summaryID >> SUMMARYID
####cat summaryM >> SUMMARYM
####cat summaryV >> SUMMARYV
cat summaryR >> SUMMARYR
cat summarySqE >> SUMMARYSqE
cat outfile >> OUTFILE
cat distributionfile >> SUMMARYq

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/data.ped $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/data.map $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/qt.phe $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/list* $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/outfileSLIM $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/data.F $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYID $WDIR/SUMMARYID_$n
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYq $WDIR/SUMMARYq_$n
####cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYM $WDIR/SUMMARYM_$n
####cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYV $WDIR/SUMMARYV_$n
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYR $WDIR/SUMMARYR_$n
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/SUMMARYSqE $WDIR/SUMMARYSqE_$n
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/OUTFILE $WDIR/OUTFILE_$n
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/checkfile $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/slimout $WDIR/
cp -r /state/partition1/slimDruet$n/$SLURM_JOBID/timefile $WDIR/

done

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/slimDruet$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
