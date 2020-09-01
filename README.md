# Caballero-et-al.-Evol-Appli-2020

Programs and software used in the simulations of 
On the estimation of inbreeding depression using different measures of inbreeding from molecular markers
Armando Caballero, Beatriz Villanueva and Tom Druet
Evolutionary Applications 2020

To compile: bash gcc
To run in an scratch directory: qsub script_SLIM_ID_100.sh <INPUT> <NIND> <REPS> <n>
INPUT = input for slim, e.g. file slimINPUT_N100
NIND = number of individuals, e.g 100
REPS = number of replicates
n = number of the scratch directory in /state/partition1/slim$n (this direction should be changed depending on the system).
