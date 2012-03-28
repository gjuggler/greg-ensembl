bsub -q research-rh6 -I -n 10 -a openmpi mpirun.lsf -np 10 -mca btl tcp,self /homes/greg/src/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < args.txt

#mpirun -np 5 /homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/bin/HYPHYMPI BASEPATH=/homes/greg/lib/greg-ensembl/projects/thesis-corrections/hyphy_build/lib/hyphy GABranchFiles/ModelSelectorBranchLocal.bf < args.txt
