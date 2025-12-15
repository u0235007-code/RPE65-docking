#! /bin/tcsh
#SBATCH --time=2:00:00 # walltime, abbreviated by -t
#SBATCH --nodes=1      # number of cluster nodes, abbreviated by -N
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm-%j.err-%N # name of the stderr, using job and first node values
#SBATCH --ntasks=1    # number of MPI tasks, abbreviated by -n
#SBATCH --account=horvath     # account - abbreviated by -A
#SBATCH --job-name=so_VFQ5
#SBATCH --partition=notchpeak-shared  # partition, abbreviated by -p
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=martin.horvath@utah.edu
#
# Update the path to the working directory which should have all the necessary input files (can be in subdirectories)
setenv WORK_DIR $HOME/Docking/dietary_lutein_2025/RPE65_dietary_lutein
setenv INPUT_DIR $WORK_DIR/INPUT_RPE65

setenv VINA_DIR $WORK_DIR/OUTPUT_RPE65_VFQ5
setenv LIGAND_DIR $HOME/../u0452756/Docking/LIGANDS/VFQ5

# Update the uID for scratchdisk space
setenv SCRATCH_DIR /scratch/general/nfs1/u0235007/myscratch$SLURM_JOBID
setenv SCOREONLY_DIR $SCRATCH_DIR/SCOREONLY_RPE65_VFQ5
setenv SPLIT_DIR $SCRATCH_DIR/SPLIT_RPE65_VFQ5
#
#
#
# Create the Scratchdisk and Output directories
mkdir -p $SCRATCH_DIR
mkdir -p $SCOREONLY_DIR
mkdir -p $SPLIT_DIR
#
# Copy files from Working directory to the Scratchdisk
cp -pR $INPUT_DIR/* $SCRATCH_DIR
cp -pR $LIGAND_DIR $SCRATCH_DIR
cp -pR $VINA_DIR $SCRATCH_DIR
#
# navigate to the Scratchdisk and run the job
cd $SCRATCH_DIR
module load autodock-vina

foreach ligand_name ( `cat $LIGAND_DIR/vfq5_names_pdbqt.csv` )
  echo -n "processing vfq5 "
  echo "$ligand_name"      
  vina_split --ligand $SPLIT_DIR/"split_RPE65_"$ligand_name"_ligand_mode_" --flex $SPLIT_DIR/"split_RPE65_"$ligand_name"_flex_mode_" --input $VINA_DIR/vina_RPE65_$ligand_name.pdbqt
  foreach mode (1 2 3 4 5)
     if (-f $SPLIT_DIR/"split_RPE65_"$ligand_name"_ligand_mode_"$mode".pdbqt") then
        vina --score_only --autobox --receptor cow_7l0e_chainA_rigid.pdbqt --ligand $SPLIT_DIR/"split_RPE65_"$ligand_name"_ligand_mode_"$mode".pdbqt" --flex $SPLIT_DIR/"split_RPE65_"$ligand_name"_flex_mode_"$mode".pdbqt" > $SCOREONLY_DIR/"split_RPE65_"$ligand_name"_mode_"$mode".log"
     endif
  end 
end 
#
# delete previous results and copy updated results to the working directory
# rm -rf $WORK_DIR/SCOREONLY_RPE65_VFQ5
# rm -rf $WORK_DIR/SPLIT_RPE65_VFQ5
cp -pR $SCOREONLY_DIR $WORK_DIR/
cp -pR $SPLIT_DIR $WORK_DIR/
#
# remove files from Scratchdisk; remove comment marker if this is safe, meaning we know the results will be copied over with confidence
rm -rf *
