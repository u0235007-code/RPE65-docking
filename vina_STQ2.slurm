#! /bin/tcsh
#SBATCH --time=8:00:00 # walltime, abbreviated by -t
#SBATCH --nodes=1      # number of cluster nodes, abbreviated by -N
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm-%j.err-%N # name of the stderr, using job and first node values
#SBATCH --ntasks=1    # number of MPI tasks, abbreviated by -n
#SBATCH --account=horvath     # account - abbreviated by -A
#SBATCH --job-name=av_STQ2
#SBATCH --partition=notchpeak  # partition, abbreviated by -p
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=martin.horvath@utah.edu
#
# Update the path to the working directory which should have all the necessary input files (can be in subdirectories)
setenv WORK_DIR $HOME/Docking/dietary_lutein_2025/RPE65_dietary_lutein
setenv LIGAND_DIR $HOME/../u0452756/Docking/LIGANDS/STQ2

setenv INPUT_DIR $WORK_DIR/INPUT_RPE65
#
# Update the uID for scratchdisk space
setenv SCRATCH_DIR /scratch/general/nfs1/u0235007/myscratch$SLURM_JOBID
#
# Update the name of the output directory which will be a subdirectory in the Working directory at the end of this job
setenv OUTPUT_DIR $SCRATCH_DIR/OUTPUT_RPE65_STQ2
#
# Create the Scratchdisk and Output directories
mkdir -p $SCRATCH_DIR
mkdir -p $OUTPUT_DIR
#
# Copy files from Working directory to the Scratchdisk
cp -pR $INPUT_DIR/* $SCRATCH_DIR
cp -pR $LIGAND_DIR $SCRATCH_DIR
#
# navigate to the Scratchdisk and run the job
cd $SCRATCH_DIR
module load autodock-vina
foreach ligand_name ( `cat STQ2/stq2_names_pdbqt.csv`)
	echo -n "processing stq2 ligand "
	echo "$ligand_name"
       	vina --cpu 12 --exhaustiveness 12 --config RPE65_config.txt --ligand STQ2/$ligand_name.pdbqt --out $OUTPUT_DIR/vina_RPE65_$ligand_name.pdbqt > $OUTPUT_DIR/vina_RPE65_$ligand_name.log
end
#
# copy results to the working directory
cp -pR $OUTPUT_DIR $WORK_DIR/
#
# remove files from Scratchdisk; remove comment marker if this is safe, meaning we know the results will be copied over with confidence
  rm -rf *
