/* This file defines the transition matrix for the Muse-Gaut 94 model x an arbitrary 4x4 rate matrix
   for nucleotide substituions.
   
   The file should be used as follows:
   
   1) Read Data File and create datafilter filteredData
   2) #include this file (or use SelectTemplateModel(filteredData);)
   3) Define the tree
   4) Proceed with the likelihood function using 'vectorOfFrequencies' as the vector to pass to the constructor.
   
   This model has the following signature:
   	#Short:MG94custom#
   	#Desc:Muse-Gaut 94 x an arbitrary 4x4 rate matrix and 9 (3x4) frequency parameters. Possible Gamma Variation.#
   	#Dimension:*#
    #DataType:codon#
   	#FileName:MG94custom.mdl#
   
   04/18/2002  by Sergei L. Kosakovsky Pond
*/

ModelMatrixDimension = 0;

global		omega = 1;
global     AC = 1;
global 	   AT = 1;
global     CG = 1;
global	   CT = 1;
global     GT = 1;



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64;
		for (h = 0 ;h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10)
			{
				ModelMatrixDimension = ModelMatrixDimension-1;
			}
		}
	}
	
	omega  = omega;
	AT = AT;
	CG = CG;
	CT = CT;
	GT = GT;
	AC = AC;

	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 
	

	hshift = 0;

	if (modelType < 2)
	{
		for (h=0; h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10) 
			{
				hshift = hshift+1;
				continue; 
			}
			vshift = hshift;
			for (v = h+1; v<64; v=v+1)
			{
				diff = v-h;
				if (_Genetic_Code[v]==10) 
				{
					vshift = vshift+1;
					continue; 
				}
				nucPosInCodon = 2;
				if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
				{
					if (h$4==v$4)
					{
						transition = v%4;
						transition2= h%4;
					}
					else
					{
						if(diff%16==0)
						{
							transition = v$16;
							transition2= h$16;
							nucPosInCodon = 0;
						}
						else
						{
							transition = v%16$4;
							transition2= h%16$4;
							nucPosInCodon = 1;
						}
					}
					if (transition<transition2)
					{
						trSM = transition;
						trLG = transition2;
					}
					else
					{
						trSM = transition2;
						trLG = transition;
					}
					
					if (trSM==0)
					{
						if (trLG==1)
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := AC*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								if (modelType)
								{
									ModelMatrixName[h-hshift][v-vshift] := AC*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AC*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := AC*nonSynRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AC*nonSynRate*EFV__[transition2__][nucPosInCodon__];
								}
							}
						}
						else
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									if (modelType)
									{
										ModelMatrixName[h-hshift][v-vshift] := omega*synRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := omega*synRate*EFV__[transition2__][nucPosInCodon__];
									}
									else
									{
										ModelMatrixName[h-hshift][v-vshift] := nonSynRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := nonSynRate*EFV__[transition2__][nucPosInCodon__];									
									}
								}							
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AT*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									if (modelType)
									{
										ModelMatrixName[h-hshift][v-vshift] := AT*omega*synRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := AT*omega*synRate*EFV__[transition2__][nucPosInCodon__];
									}
									else
									{
										ModelMatrixName[h-hshift][v-vshift] := AT*nonSynRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := AT*nonSynRate*EFV__[transition2__][nucPosInCodon__];									
									}
								}							
							}
						}
					}
					else
					{
						if (trSM==1)
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CG*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									if (modelType)
									{
										ModelMatrixName[h-hshift][v-vshift] := CG*omega*synRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := CG*omega*synRate*EFV__[transition2__][nucPosInCodon__];
									}
									else
									{
										ModelMatrixName[h-hshift][v-vshift] := CG*nonSynRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := CG*nonSynRate*EFV__[transition2__][nucPosInCodon__];									
									}
								}
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CT*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									if (modelType)
									{
										ModelMatrixName[h-hshift][v-vshift] := CT*omega*synRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := CT*omega*synRate*EFV__[transition2__][nucPosInCodon__];
									}
									else
									{
										ModelMatrixName[h-hshift][v-vshift] := CT*nonSynRate*EFV__[transition__][nucPosInCodon__];
										ModelMatrixName[v-vshift][h-hshift] := CT*nonSynRate*EFV__[transition2__][nucPosInCodon__];
									}
									
								}							
							}
						}
						else
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := GT*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								if (modelType)
								{
									ModelMatrixName[h-hshift][v-vshift] := GT*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := GT*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := GT*nonSynRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := GT*nonSynRate*EFV__[transition2__][nucPosInCodon__];								
								}
							}							
						}
					}
				}
		   }
	    }		
	}
	else
	{
		for (h=0; h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10) 
			{
				hshift = hshift+1;
				continue; 
			}
			vshift = hshift;
			for (v = h+1; v<64; v=v+1)
			{
				diff = v-h;
				if (_Genetic_Code[v]==10) 
				{
					vshift = vshift+1;
					continue; 
				}
				nucPosInCodon = 2;
				if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
				{
					if (h$4==v$4)
					{
						transition = v%4;
						transition2= h%4;
					}
					else
					{
						if(diff%16==0)
						{
							transition = v$16;
							transition2= h$16;
							nucPosInCodon = 0;
						}
						else
						{
							transition = v%16$4;
							transition2= h%16$4;
							nucPosInCodon = 1;
						}
					}
					if (transition<transition2)
					{
						trSM = transition;
						trLG = transition2;
					}
					else
					{
						trSM = transition2;
						trLG = transition;
					}
					
					if (trSM==0)
					{
						if (trLG==1)
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*c*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := AC*c*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*c*omega*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := AC*c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
							}
						}
						else
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := c*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := c*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := c*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*c*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AT*c*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*c*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AT*c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
						}
					}
					else
					{
						if (trSM==1)
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*c*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CG*c*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*c*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CG*c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*c*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CT*c*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*c*omega*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CT*c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
						}
						else
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*c*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := GT*c*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*c*omega*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := GT*c*omega*synRate*EFV__[transition2__][nucPosInCodon__];
							}							
						}
					}
				}
		   }
	    }	
	}	
	if (Abs(MGCustomModelConstraintString))
	{
		ExecuteCommands (MGCustomModelConstraintString);
	}
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

categoriesUsed = 0;

if (!SKIP_MODEL_PARAMETER_LIST)
{
	#include "modelParameters.mdl";
}

sharedFlag = 1;

if (modelType > 1)
{
	categoriesUsed = 1;
	if (modelType == 2)
	{
		#include "defineGamma.mdl";
	}
	if (modelType == 3)
	{
		#include "defineHM.mdl";
	}
}

if (!SKIP_MODEL_PARAMETER_LIST)
{
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}	
}
modelDesc = "010010";
			
MGCustomRateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
paramCount	  = 0;

MGCustomModelConstraintString = "";

for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
		{
			if (MGCustomRateBiasTerms[customLoopCounter2] == "1")
			{
				MGCustomModelConstraintString = MGCustomModelConstraintString + MGCustomRateBiasTerms[customLoopCounter]+":="+MGCustomRateBiasTerms[customLoopCounter2]+";";
			}
			else
			{
				MGCustomModelConstraintString = MGCustomModelConstraintString + MGCustomRateBiasTerms[customLoopCounter2]+":="+MGCustomRateBiasTerms[customLoopCounter]+";";			
			}
			break;
		}
	}
}				

HarvestFrequencies (observedFreq,codon_dsf,3,1,1);

NICETY_LEVEL = 3;

MG94custom = 0;

MULTIPLY_BY_FREQS = PopulateModelMatrix ("MG94custom", observedFreq);

FREQUENCY_SENSITIVE = 1;

vectorOfFrequencies = BuildCodonFrequencies (observedFreq);

Model MG94customModel = (MG94custom,vectorOfFrequencies,0);

USE_POSITION_SPECIFIC_FREQS = 1;
