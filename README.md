RPE65 was docked with different classes of ligands: LUT as found in the Protein Data Bank, dietary lutein (L), meso-zeaxanthin (MZ), dietary zeaxanthin (Z), beta-carotene, and the chemical compound difluoro-emisustat. L, MZ, and Z belong to the family of molecules called xanthophylls. L, MZ, Z, and beta-carotene belong to the family of molecules called carotenoids. The molecular docking was accomplished with AutoDock VINA. A research artical describing the molecules and results of docking will be submitted for review DEC 2025.

Branches in this repository are for the different molecules and the docking outcomes. Programs for generating the molecules and docking outomes are in the approprate branch.

LUT, The lutein-like ligands as found in PDB structures of photosystems and light harvesting complexes with incorrect 3S stereo-configuration for the beta ionone ring. Dietary (3R, 3'R, 6'R) lutein has 3R configuration for the beta ionone ring. Non-dietary (3S, 3'R, 6'R) LUT ligands in the PDB have 3S configuration for the beta ionone ring.

LUT/LUT2, Lutein-like ligands with incorrect 3S beta stereo-configuration, with charges and atom types added by AutoDock Tools and no rotating bonds (rigid ligands).

LUT/LUT5, Lutein-like ligands with incorrect 3S beta stereo-configuration, with charges and atom types added by AutoDock Tools and rotating bonds connecting the ionone rings to the polyene chain (two rotating bonds).

BCR, Beta carotene molecules as found in the PDB structures of photosystems and light harvesting complexes.

BCR/BCR2, Beta carotene ligands with charges and atom types added by AutoDock Tools and no rotating bonds (rigid ligands).

BCR/BCR5, Beta carotene ligands with charges and atom types added by AutoDock Tools and rotating bonds connecting the ionone rings to the polyene chain (two rotating bonds).

STQ, Dietary (3R, 3'R, 6'R) lutein derived from lutein-like LUT ligands by reflecting the positions of atoms in the beta ionone ring across the plane defined by the double bond connecting atoms C6 and C5. The resulting (3R, 3'R, 6'R) lutein molecules are the molecules common to the human diet. The vector math for the reflection operation needed to convert 3S beta to 3R beta is encoded in the R program dietary.lutein.R and is summarized here

Briefly (x,y,z) coordinates for each atom in the beta ionone ring (C1, C2, C3, O3, C4, C5, C6, C16, C17, C18) defined corresponding vectors. The cross product (C1 - C6) x (C5 - C6) defined R, a vector perpendicular to the reflecting plane defined by atoms C1, C5 and C6 (Equation S1). 

	        R = (C1 - C6) x (C5 - C6)	Eq S1
            
For each atom being reflected, its original position defined a difference vector relative to an atom in the reflecting plane, either C1 for the two methyl groups C17 and C18, or C6 for all other atoms. The dot product of this difference vector with R provided the distance from the reflecting plane, and thereby a means to compute the projection vector along R. The coordinates of each atom could then be transformed by subtracting this projection vector scaled by 2. Equations below illustrate these computations for reflection of atom O3. Similar calculations were applied to all atoms in the beta ionone ring, except those atoms that are already in the reflecting plane (C1, C4, C5, and C18). 

	        distance(O3)  = R • (O3 - C6)	Eq. S2
	        projection(O3) = R distance(O3) / R • R	Eq. S3
	        reflection(O3) = O3 - 2 projection(O3)	Eq. S4

STQ/STQ2, Dietary (3R, 3'R, 6'R) lutein derived from lutein-like LUT ligands by reflecting the positions of atoms in the beta ionone ring across the plane defined by the double bond connecting atoms C6 and C5, with charges and atom types added by AutoDock Tools and no rotating bonds (rigid ligands).

STQ/STQ5, Dietary (3R, 3'R, 6'R) lutein derived from lutein-like LUT ligands by reflecting the positions of atoms in the beta ionone ring across the plane defined by the double bond connecting atoms C6 and C5, with charges and atom types added by AutoDock Tools and rotating bonds connecting the ionone rings to the polyene chain (two rotating bonds).

VFQ, (3R, 3'S meso-) Zeaxanthin derived from dietary lutein and lutein-like LUT ligands by replacing the epsilon ionone ring found in STQ molecules with the 3S beta ring found in LUT molecules. The resulting (3R, 3'S meso-)zeaxanthin molecules are only rarely found in the human diet. A major question for carotenoid biochemistry is what is the source of meso-zeaxanthin found in the retinas of humans and birds. Replacement of the epsilon ionone ring with the 3S beta ionone ring was accomplished by superposition of two xanthophyll molecules, one STQ and one LUT. 

In these epsilon ring replacement operations, the second beta ring was positioned as closely as possible with the epsilon ring by superposition of C7, C8 and C9 atoms of the polyene chain of the LUT molecule onto C27, C28 and C29 atoms in the STQ molecule. Coordinates for the hybrid ligand molecule were written by selectively merging C1, C2, C3, O3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16, C17, C18, C19, C20, C29, C30, C31, C32, C33, C34, C35, C39, C40 atoms belonging to the beta ring and polyene chain from the STQ molecule, thereby excluding the atoms belonging to the epsilon ionone ring, with the C1, C2, C3, O3, C4, C5, C6, C7, C8, C16, C17, C18 atoms from the beta ionone ring of the LUT molecule after superposition. Superposition and the merging of atoms was facilitated by ChimeraX scripts, superimpose.Sbeta.to.epsilon.cxc for MZ ligands. Further cleanup of atom names was accomplished with a tcsh script cleanup.mesozeaxanthin.script and a set of vi commands encoded with cleanup.mesozeaxanthin.vi for the MZ ligands. These scripts for superposition and for atom name cleanup were created by an R program dietary.zeaxanthin.R 

VFQ/VFQ2, (3R, 3'S meso-) Zeaxanthin derived from dietary lutein and lutein-like LUT ligands by replacing the epsilon ionone ring found in STQ molecules with the 3S beta ring found in LUT molecules, with charges and atom types added by AutoDock Tools and no rotating bonds (rigid ligands).

VFQ/VFQ5, (3R, 3'S meso-) Zeaxanthin derived from dietary lutein and lutein-like LUT ligands by replacing the epsilon ionone ring found in STQ molecules with the 3S beta ring found in LUT molecules, with charges and atom types added by AutoDock Tools and rotating bonds connecting the ionone rings to the polyene chain (two rotating bonds).

8UH, Dietary (3R, 3'R) zeaxanthin derived by merging two dietary lutein molecules to replacing the epsilon ionone ring found in one STQ molecules with the 3R beta ring found in a second STQ molecule. The resulting (3R, 3'R) zeaxanthin molecules are found commonly in the human diet.

In these epsilon ring replacement operations, the second beta ring was positioned as closely as possible with the epsilon ring by superposition of C7, C8 and C9 atoms of the polyene chain of the second STQ molecule onto C27, C28 and C29 atoms in the first STQ molecule. Coordinates for the hybrid ligand molecule were written by selectively merging C1, C2, C3, O3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16, C17, C18, C19, C20, C29, C30, C31, C32, C33, C34, C35, C39, C40 atoms belonging to the beta ring and polyene chain from the first STQ molecule, thereby excluding the atoms belonging to the epsilon ionone ring, with the C1, C2, C3, O3, C4, C5, C6, C7, C8, C16, C17, C18 atoms from the beta ionone ring of the second STQ molecule after superposition. Superposition and the merging of atoms was facilitated by ChimeraX scripts superimpose.beta.to.epsilon.cxc for Z ligands. Further cleanup of atom names was accomplished with a tcsh script cleanup.zeaxanthin.script and a set of vi commands encoded with cleanup.zeaxanthin.vi for the Z ligands. These scripts for superposition and for atom name cleanup were created by an R program dietary.zeaxanthin.R 

8UH/8UH2, Dietary (3R, 3'R) zeaxanthin derived by merging two dietary lutein molecules to replacing the epsilon ionone ring found in one STQ molecules with the 3R beta ring found in a second STQ molecule, with charges and atom types added by AutoDock Tools and no rotating bonds (rigid ligands).

8UH/8UH5, Dietary (3R, 3'R) zeaxanthin derived by merging two dietary lutein molecules to replacing the epsilon ionone ring found in one STQ molecules with the 3R beta ring found in a second STQ molecule, with charges and atom types added by AutoDock Tools and rotating bonds connecting the ionone rings to the polyene chain (two rotating bonds).


