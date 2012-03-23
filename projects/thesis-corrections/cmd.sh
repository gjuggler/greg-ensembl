bsub -q research-rh6 -Is -n 5 -a openmpi "mpirun.lsf -mca btl tcp,self /homes/greg/src/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < args.txt"

#mpirun -np 5 /homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < args.txt

#mpirun -np 1 /homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/BranchGAResultProcessor.bf < args2.txt