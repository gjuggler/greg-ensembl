Model_Name  		= "GY94_3x4";
Model_Options 		= 3;
Model_Dimension 	= 64;
Model_EFV_Type		= "Observed Codon";

ModelMatrixDimension = 0;

rateClassCount = 3;

function BuildCodonFrequencies (EFV)
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
			PIStop = PIStop-EFV[first][0]*EFV[second][1]*EFV[third][2];
			continue; 
		}
		result[h-hshift][0]=EFV[first][0]*EFV[second][1]*EFV[third][2];
	}
	return result*(1.0/PIStop);
}

function PopulateModelMatrix (ModelMatrixName&, EFV)
{	
	global globalVariable_TVTS = 0.1;

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
	
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 

	hshift = 0;

	if (modelType == 0)
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
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
			  		{
			  			if ((Abs(transition-transition2)%2)==1)
			  			{
			  				ModelMatrixName[h-hshift][v-vshift] := globalVariable_TVTS*synRate;
			  				ModelMatrixName[v-vshift][h-hshift] := globalVariable_TVTS*synRate;
			  			
			  			}
			  			else
			  			{
			  				ModelMatrixName[h-hshift][v-vshift] := synRate;
			  				ModelMatrixName[v-vshift][h-hshift] := synRate;
			  			}
				  	}
			  		else
			  		{
			  			if ((Abs(transition-transition2)%2)==1)
			  			{
			  				ModelMatrixName[h-hshift][v-vshift] := globalVariable_TVTS*nonSynRate;
			  				ModelMatrixName[v-vshift][h-hshift] := globalVariable_TVTS*nonSynRate;
			  			
			  			}
			  			else
			  			{
			  				ModelMatrixName[h-hshift][v-vshift] := nonSynRate;
			  				ModelMatrixName[v-vshift][h-hshift] := nonSynRate;
			  			}
				  	}
			  	}
			}
		}
	}
	else
	{
		global omega = 1;

		if (modelType == 1)
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
				  			}
				  			else
				  			{
				  				transition = v%16$4;
				  				transition2= h%16$4;
				  			}
				  		}
				  		if (_Genetic_Code[h]==_Genetic_Code[v]) 
				  		{
				  			if ((Abs(transition-transition2)%2)==0)
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := synRate;
				  			
				  			}
				  			else
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := globalVariable_TVTS*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := globalVariable_TVTS*synRate;
				  			}
					  	}
				  		else
				  		{
				  			if ((Abs(transition-transition2)%2)==0)
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := omega*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := omega*synRate;
				  			
				  			}
				  			else
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := omega*globalVariable_TVTS*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := omega*globalVariable_TVTS*synRate;
				  			}
					  	}
				  	}
			    }
			}	
		}
		else
		{
			global shapeParameter = .5;
			shapeParameter:>0.01;
			shapeParameter:<100;
			category     categoryVariable = 
						(rateClassCount, EQUAL, MEAN, GammaDist(_x_,shapeParameter,shapeParameter), CGammaDist(_x_,shapeParameter,shapeParameter), 0 , 
				  									  1e25,CGammaDist(_x_,shapeParameter+1,shapeParameter));
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
				  			}
				  			else
				  			{
				  				transition = v%16$4;
				  				transition2= h%16$4;
				  			}
				  		}
				  		
				  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
				  		{
				  			if ((Abs(transition-transition2)%2)==1)
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := categoryVariable*globalVariable_TVTS*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := categoryVariable*globalVariable_TVTS*synRate;
				  			
				  			}
				  			else
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := categoryVariable*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := categoryVariable*synRate;
				  			}
					  	}
				  		else
				  		{
				  			if ((Abs(transition-transition2)%2)==0)
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := categoryVariable*omega*globalVariable_TVTS*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := categoryVariable*omega*globalVariable_TVTS*synRate;
				  			
				  			}
				  			else
				  			{
				  				ModelMatrixName[h-hshift][v-vshift] := categoryVariable*omega*synRate;
				  				ModelMatrixName[v-vshift][h-hshift] := categoryVariable*omega*synRate;
				  			}
					  	}
				  	}
				 }
			}		
		}
	}
	
	return 1;
}


if (!SKIP_MODEL_PARAMETER_LIST)
{
	#include "modelParameters.mdl";
}
if (modelType > 1)
{
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

HarvestFrequencies (observedFreq,codon_dsf,3,1,1);

NICETY_LEVEL = 3;

GY94custom = 0;

MULTIPLY_BY_FREQS = PopulateModelMatrix ("GY94custom", observedFreq);

FREQUENCY_SENSITIVE = 1;

vectorOfFrequencies = BuildCodonFrequencies (observedFreq);

Model GY94customModel = (GY94custom,vectorOfFrequencies,0);

USE_POSITION_SPECIFIC_FREQS = 1;
