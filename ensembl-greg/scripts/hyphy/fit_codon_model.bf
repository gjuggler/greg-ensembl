/* See http://www.hyphy.org/pubs/Nef/galaxy/ */

NICETY_LEVEL = 3;
VERBOSITY_LEVEL = 0;

skipCodeSelectionStep = 1;
#include "chooseGeneticCode.def";

DataSet codon_ds = ReadDataFile ("test.fa");
DataSetFilter codon_dsf = CreateFilter (codon_ds,3,"","",GeneticCodeExclusions);
fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n", codon_ds,"\n");

/* Set model parameters here. */

SKIP_MODEL_PARAMETER_LIST = 1;
modelType = 1;
SKIP_DISTRIBUTION_CHOICE = 1;
distribution = 0;
SKIP_DISTRIBUTION_CLASSES = 1;
resp = 2;

/* generate codon model (MG94customModel) */
#include "fit_codon_model.ibf";

/* read in tree from file */
ACCEPT_BRANCH_LENGTHS = 1;
ACCEPT_ROOTED_TREES = 1;
fscanf ("test.nh", "Raw", tree_string);
fprintf(stdout,"Tree: ",tree_string,"\n");
Tree	codon_tree = tree_string;


LikelihoodFunction codon_lf = (codon_dsf, codon_tree);

AUTO_PARALLELIZE_OPTIMIZE = 1;	/* attempt to use MPI */
Optimize (res, codon_lf);
AUTO_PARALLELIZE_OPTIMIZE = 0;

COVARIANCE_PARAMETER = "R";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx,codon_lf);
fprintf (stdout, "omega_lo:", cmx[0], "\n");
fprintf (stdout, "omega:", cmx[1], "\n");
fprintf (stdout, "omega_hi:", cmx[2], "\n");

LIKELIHOOD_FUNCTION_OUTPUT = 7;		/* save LF to file */
fprintf("test.out", codon_lf);
LIKELIHOOD_FUNCTION_OUTPUT = 2;		/* reset to default (?) */

fprintf (stdout, "\n______________RESULTS______________\n",codon_lf);



/* UNUSED CRAP: */
/*
function computeScalingFactorB(rateMatrix, baseFreqs)
{
	B = 0;
	for (n1 = 0; n1 < Rows(rateMatrix); n1 = n1+1)
	{
		for (n2 = 0; n2 < Columns(rateMatrix); n2 = n2+1)
		{
			if (n2 != n1)
			{
				B = B + baseFreqs[n1]*baseFreqs[n2]*rateMatrix[n1][n2];
			}
		}
	}
	return B;
}
global scalingB	= computeScalingFactorB (MG94custom, vectorOfFrequencies);
branchNames 	= BranchName (codon_tree, -1);
branchLengths	= BranchLength (codon_tree, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("codon_tree." + branchNames[k] + ".synRate:=" + branchLengths[k] + "/scalingB;");
}
*/
