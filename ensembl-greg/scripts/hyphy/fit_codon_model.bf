/* See http://www.hyphy.org/pubs/Nef/galaxy/ */
NICETY_LEVEL = 3;
VERBOSITY_LEVEL = 0;

skipCodeSelectionStep = 1;
#include "chooseGeneticCode.def";

SKIP_MODEL_PARAMETER_LIST = 1;
SKIP_DISTRIBUTION_CHOICE = 1;
SKIP_DISTRIBUTION_CLASSES = 1;

/* Give us two filenames on standard input. */
fscanf (stdin, "String",  alignmentFilename);
fscanf (stdin, "String",  treeFilename);
fscanf (stdin, "String",  modelFilename);

fprintf(stdout,"\n",alignmentFilename,"\n",treeFilename,"\n",modelFilename,"\n");
DataSet codon_ds = ReadDataFile (alignmentFilename);
DataSetFilter codon_dsf = CreateFilter (codon_ds,3,"","",GeneticCodeExclusions);
fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n", codon_ds,"\n");

/* Set model parameters here. */
/*
Possible modelType values:
 0: Local (per-branch optimized)
 1: Global (global throughout tree -- good default)
 2: Global w/variation
 3: Global w/variation + HM
*/
modelType = 1;
/*
If modelType is >= 2, we need a distribution:
 0:Gamma, 1:Beta, 2:Normal(x>0), 3:Lognormal, 4:Gamma+Inv, 5:GammaQuadratures, 6:InvClass, 7:GeneralDiscrete,
 8:BetaGamma, 9:BetaLognormal, 10:Normal(x>0)ModBeta, 11:BetaUniform, 12:Uniform, 13:BetaGammaInv
*/
distribution = 3;
/* Number of rate class bins */
resp = 3;

/* generate codon model */
/*#include "gy94.bf";*/
#include "mg94.bf";

/* read in tree from file */
ACCEPT_BRANCH_LENGTHS = 1;
ACCEPT_ROOTED_TREES = 1;
fscanf (treeFilename, "Raw", tree_string);
fprintf(stdout,"Tree: ",tree_string,"\n");
Tree	codon_tree = tree_string;

LikelihoodFunction codon_lf = (codon_dsf, codon_tree);

AUTO_PARALLELIZE_OPTIMIZE = 1;	/* attempt to use MPI */
Optimize (res, codon_lf);
AUTO_PARALLELIZE_OPTIMIZE = 0;

COVARIANCE_PARAMETER = "omega";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx,codon_lf);
fprintf (stdout, "omega_lo:", cmx[0], "\n");
fprintf (stdout, "omega:", cmx[1], "\n");
fprintf (stdout, "omega_hi:", cmx[2], "\n");

LIKELIHOOD_FUNCTION_OUTPUT = 7;		/* save LF to file */
/*fprintf("test.out", codon_lf);*/
LIKELIHOOD_FUNCTION_OUTPUT = 2;		/* reset to default (?) */

fprintf (stdout, "\n______________RESULTS______________\n",codon_lf);
