#include <R.h>
#include <math.h>
#include <Rmath.h>

// #define WITH_PERMUTATIONS
#define FAST_VERSION_INFOREP

// #define DEBUG_TOOL
#ifdef DEBUG_TOOL
	FILE *file;

	#define PRINT(_type, _x) \
		{ \
			fprintf(file, #_x" "_type"\r\n", _x); \
			fflush(file); \
		}
		
	#define PRINT_SET_FUNCTION(_type, _set, _lenSet) \
		{ \
			fprintf(file, #_set" "); \
			int i = 0; \
			for(i = 0; i < (_lenSet); ++i) \
			{ \
				fprintf(file, _type" ", _set[i]); \
			} \
			fprintf(file, "\r\n"); \
			fflush(file); \
		}
#endif

#define MALLOC(_type, _size) (_type*)malloc((_size)*sizeof(_type))


int Powi(int x, int y)
{
	int tmp = 1;
	while(y)
	{
		tmp = tmp*x;
		y--;
	}
	return tmp;
}

int Fact(int n)
{
	int tmp = 1;
	while(n > 1)
	{
		tmp = tmp*n;
		--n;
	}
	return tmp;
}

int Contains(double _elt, double *_set, int _setSize)
{
	for(int i = 0; i < _setSize; ++i)
	{
		if(_set[i] == _elt)
		{
			return 1;
		}
	}
	return 0;
}


double FindFrequence(double value, double *alleleDouble, double *frequence, int lenFreq)
{
	int i;
	for(i = 0; i < lenFreq; ++i)
	{
		if(alleleDouble[i] == value)
		{
			return frequence[i];
		}
	}
	return 1.;
}


void infoRepC3(double *R, int *lenR, double *Ggi, int *lenGgi, double *prDHet, int *lenHet, double *prDHom, int *lenHom, double *prC, double *q, double *frequence, int *lenFreq, double *res)
{
	#ifdef DEBUG_TOOL
		//file = fopen("debug2.txt","w+b");
		fprintf(file, "begin debug2\r\n");
		fflush(file);
		if(R != NULL && *lenR >= 3)
		{
			fprintf(file, "res = %f, R=%f\r\n", *res, *R);
			fflush(file);
			fprintf(file, "R1 = %f, R3 =%f\r\n", R[0], R[2]);
			fflush(file);
		}
		// fprintf(file, "allele = %s, allele=%s\r\n", allele[0], allele[1]);
		// fflush(file);
		fprintf(file, "frequence = %f, %f\r\n", frequence[0], frequence[1]);
		fflush(file);
	#endif

	int i = 0;
	int j = 0;
	int k = 0;
	//int l = 0;
	
	double *alleleDouble = q;
	
	// create union set between R and Ggi
	double *UnionRAndGgi = MALLOC(double,*lenR + *lenGgi);
	int lenUnionRAndGgi = 0;
	
	if(R != NULL)
	{
		memcpy(UnionRAndGgi, R, *lenR*sizeof(double));
		lenUnionRAndGgi = *lenR;
	}
	
	for(i = 0; i < *lenGgi; ++i)
	{
		for(j = 0; j < lenUnionRAndGgi; ++j)
		{
			if(UnionRAndGgi[j] == Ggi[i])
			{
				break;
			}
		}
		if(j == lenUnionRAndGgi)
		{
			UnionRAndGgi[lenUnionRAndGgi] = Ggi[i];
#ifdef DEBUG_TOOL
			fprintf(file, "UnionRAndGgi[%d] %f\r\n", lenUnionRAndGgi, UnionRAndGgi[lenUnionRAndGgi]);
			fflush(file);
#endif
			lenUnionRAndGgi++;
		}
	}
#ifdef DEBUG_TOOL
	fprintf(file, "lenUnionRAndGgi %d\r\n", lenUnionRAndGgi);
	fflush(file);
#endif

	// find good probabilities to use for each allele
	double *prD = MALLOC(double, *lenGgi/2);
	for(i = 0; i < *lenGgi; i += 2)
	{
		prD[i/2] = prDHet[i/2];
		// if homozygote
		if(Ggi[i] == Ggi[i + 1])
		{
			prD[i/2] = prDHom[i/2];
#ifdef DEBUG_TOOL
			fprintf(file, "prD[%d] %f\r\n", i/2, prD[i/2]);
			fflush(file);
#endif
		}
	}
	
	double product = 1.;
	int isContaminant = 0;
	for(i = 0; i < lenUnionRAndGgi; ++i)
	{
#ifdef DEBUG_TOOL
		fprintf(file, "for each x in unionRandGgi %d %f\r\n", i, UnionRAndGgi[i]);
		fflush(file);
#endif
		double sum = 0.;
		
		double *proba = MALLOC(double, (*lenGgi)/2);
		int lenProba = 0;
		for(j = 0; j < (*lenGgi); j += 2)
		{
			if(UnionRAndGgi[i] == Ggi[j] || UnionRAndGgi[i] == Ggi[j + 1])
			{
				proba[lenProba] = prD[j/2];
#ifdef DEBUG_TOOL
				fprintf(file, "proba[%d] %f\r\n", lenProba, proba[lenProba]);
				fflush(file);
#endif
				lenProba++;
			}
		}
		
		// if allele is in R
		if(R != NULL && Contains(UnionRAndGgi[i], R, *lenR) == 1)
		{
			// if allele is not a contaminant
			if(lenProba != 0)
			{
#ifdef FAST_VERSION_INFOREP
				double p = 1.;
				for(k = 0; k < lenProba; ++k)
				{
					p *= proba[k];
				}
				sum += 1. - p;
#else
				int pow = Powi(2, lenProba);
				for(j = 1; j < pow; ++j)
				{
					double p = 1.;
					for(k = 0; k < lenProba; ++k)
					{
						double drop = proba[k];
						if(j & (1 << k))
						{
							drop = 1. - drop;
						}
						p *= drop;
					}
					sum += p;
#ifdef DEBUG_TOOL
					fprintf(file, "sum += p %f\r\n", sum);
					fflush(file);
#endif
				}
#endif
			}
			else
			{
				sum = (*prC)*FindFrequence(UnionRAndGgi[i], alleleDouble, frequence, *lenFreq);
				isContaminant = 1;
#ifdef DEBUG_TOOL
				fprintf(file, "isContaminant sum = %f\r\n", sum);
				fflush(file);
#endif
			}
		}
		// if allele is not in R
		else
		{
			double p = 1.;
			for(j = 0; j < lenProba; ++j)
			{
				p *= proba[j];
			}
			sum = p;
#ifdef DEBUG_TOOL
			fprintf(file, "sum = p %f\r\n", sum);
			fflush(file);
#endif
		}
		free(proba);
		product *= sum;
	}
	
	if(isContaminant == 0)
	{
		product *= (1. - *prC);
	}
	
	*res = product;
	
#ifdef DEBUG_TOOL
    fprintf(file, "end debug2\r\n");
#endif

	free(prD);
	free(UnionRAndGgi);
}


double infoRepLoop2(double *_repList, int _lenRepList, double *_ggi, int _lenGgi, double *_prDHet, int _lenHet
				, double *_prDHom, int _lenHom, double _prC, double *_q, double *_frequence, int _lenFreq)
{
	double prodResRep = 1.;
	int inc = 0;
	while(inc < _lenRepList)
	{
		double *CurRepliste = &(_repList[inc]);
		int lenCurRepliste = 0;
		while(lenCurRepliste < (_lenRepList - inc) && _repList[inc + lenCurRepliste] != 0.)
		{
			lenCurRepliste = lenCurRepliste + 1;
		}

#ifdef DEBUG_TOOL
		PRINT_SET_FUNCTION("%f", CurRepliste, lenCurRepliste);
		PRINT_SET_FUNCTION("%f", _ggi, _lenGgi);
#endif
		
		double resRep;
		//void infoRepC3(double *R, int *lenR, double *Ggi, int *lenGgi, double *prDHet, int *lenHet, double *prDHom, int *lenHom, double *prC, double *q, double *frequence, int *lenFreq, double *res)
		infoRepC3(CurRepliste, &lenCurRepliste, _ggi, &_lenGgi, _prDHet, &_lenHet, _prDHom, &_lenHom
							, &_prC, _q, _frequence, &_lenFreq, &resRep);
							
#ifdef DEBUG_TOOL
		PRINT("%f", resRep);
#endif

		prodResRep = prodResRep*resRep;
		
		inc = inc + lenCurRepliste + 1;
		if(lenCurRepliste == 0)
		{
			inc++;
		}
	}
	
	return prodResRep;
}


//function to calcule the probability of Pr(T,U,V|H); but used as: Pr(U|T,V,H)
double Pevid8(double *_uset, int _lenUset, double *_q, double *_frequence, int _lenFreq
			, double *_t, int _lenT, double *_v, int _lenV, double _theta)
{
	if(_lenT == 1)
	{
		_lenT = 0;
	}
	
	if(_lenV == 1)
	{
		_lenV = 0;
	}
	
	//file = fopen("debug3.txt","w+b");
	// fprintf(file, "begin debug3\r\n");
	// fflush(file);
	double *TuV = MALLOC(double, _lenT + _lenV + _lenUset);
	int lenTuV = _lenT + _lenV;
	memcpy(TuV, _t, _lenT*sizeof(double));
	memcpy(&(TuV[_lenT]), _v, _lenV*sizeof(double));
	
	int i;
	int j;
	double prod = 1;
	
	for(i = 0; i < _lenUset; ++i)
	{
		double index = 1.;
		if(i%2 == 0)
		{
			if(_uset[i] == _uset[i + 1]) // is homozygote _frequence[_uset[i]]^2
			{
				//index = 1.;
			}
			else // *_frequence[_uset[i]]*_frequence[_uset[i+1]])
			{
				index = 2.;
			}
		}
		
		// compute how many time we have _uset[i] in Tuv
		// number of already observed alleles of type m
		int m = 0;
		for(j = 0; j < lenTuV; ++j)
		{
			if(_uset[i] == TuV[j])
			{
				m++;
			}
		}
		int n = lenTuV; // number of already sampled alleles
		
		// find frequence
		double freq = 0.;
		for(j = 0; j < _lenFreq; ++j)
		{
			if(_q[j] == _uset[i])
			{
				freq = _frequence[j];
				break;
			}
		}
		
		double prop = (m*_theta + (1 - _theta)*freq)/(1 + (n - 1)*_theta) ;
		
		TuV[lenTuV] = _uset[i];
		lenTuV++;
		
		prod = prod*prop*index;
	}
	
	free(TuV);

#ifdef WITH_PERMUTATIONS
	return prod;
#else
	// compute nb of different permutations pairs	
	int *nbPairs = MALLOC(int, (_lenUset/2 + 1));	
	/*for(i = 0; i < _lenUset/2 + 1; ++i)
	{
		nbPairs[i] = 0;
	}*/
	memset(nbPairs, 0, sizeof(int)*(_lenUset/2 + 1));
	double *diffPairs = MALLOC(double, _lenUset);	
	for(i = 0; i < _lenUset; i = i + 2)
	{
		for(j = 0; nbPairs[j] != 0; ++j)
		{
			// WARNING it works if we are sure that a pair can be only [20,21] and never [21,20]
			if(_uset[i] == diffPairs[j*2] && _uset[i + 1] == diffPairs[j*2 + 1])
			{
				break;
			}
		}
		nbPairs[j] = nbPairs[j] + 1;
		diffPairs[j*2] = _uset[i];
		diffPairs[j*2 + 1] = _uset[i + 1];

#ifdef DEBUG_TOOL
		fprintf(file, "diffPairs[j*2] %d\r\n", diffPairs[j*2]);
#endif
	}
	
	int prodPermutation = 1;
	int nb = 0;
	for(j = 0; nbPairs[j] != 0; ++j)
	{
#ifdef DEBUG_TOOL
		fprintf(file, "nbPairs[j] %d\r\n", nbPairs[j]);
		fflush(file);
#endif

		prodPermutation = prodPermutation*Fact(nbPairs[j]);
		nb = nb + nbPairs[j];
	}
	int nbPermutations = Fact(nb)/prodPermutation;

	free(nbPairs);
	free(diffPairs);
	
#ifdef DEBUG_TOOL
		fprintf(file, "nb %d\r\n", nb);
		fprintf(file, "nbPermutations %d\r\n", nbPermutations );
		fflush(file);
#endif

	return nbPermutations*prod;
#endif
}

// double *Repliste, int *lenRepliste: vector of alleles observed at a given replicate,  different replicates are seperated by a 0
// double *T, int *lenT: genotypes of the known contributors, given as a vector of size [2]
// double *V, int *lenV: genotypes of the known non-contributors, given as a vector of size [2]
// int *x: number of unknown contributors
// double *theta: scalare, theta correction
// double *prDHet, int *lenHet: dropout of heterozygotes, vector, must of size (T + x*2)
// double *prDHom, int *lenHom: dropout of homozygote, vector, must of size (T + x*2)
// double *prC: scalar
// char **allele, double *frequence, int *lenFreq: allele is a char given the list of alleles at a given locus, with frequencies *frequence
// double *sortieR: output for R 
void evidenceC3(double *Repliste, int *lenRepliste,double *T, int *lenT,double *V, int *lenV, int *x, double *theta,
double *prDHet, int *lenHet, double *prDHom, int *lenHom, double *prC, char **allele, double *frequence, int *lenFreq, double *sortieR)
{
#ifdef DEBUG_TOOL
	file = fopen("debug.txt", "w+b");
	fprintf(file, "begin debug\r\n");
	fflush(file);
	PRINT_SET_FUNCTION("%f", Repliste, *lenRepliste);
	PRINT_SET_FUNCTION("%f", T, *lenT);
	PRINT_SET_FUNCTION("%f", V, *lenV);
	PRINT_SET_FUNCTION("%f", prDHet, *lenHet);
	PRINT_SET_FUNCTION("%f", prDHom, *lenHom);
	PRINT_SET_FUNCTION("%s", allele, *lenFreq);
	PRINT_SET_FUNCTION("%f", frequence, *lenFreq);
#endif

	int i;
	int j;
	
	double *Q = MALLOC(double, *lenFreq);
	for(i = 0; i < *lenFreq; ++i)
	{
		Q[i] = atof(allele[i]);
	}
		
	// we compute combination only if x different of 0
	if(*x != 0)
	{
		// combinaison avec remise, definit les genotypes possible pour une seule personne
		// create unique pair set
		int nbAllelePair = 0;
		for(i = 1; i <= *lenFreq; ++i)
		{
			nbAllelePair = nbAllelePair + i;
		}

#ifdef DEBUG_TOOL
		fprintf(file, "nbAllelePair %d\r\n", nbAllelePair);
		fflush(file);
#endif
		
		int *allelePair = MALLOC(int, nbAllelePair*2);
		int inc = 0;
		for(i = 0; i < *lenFreq; ++i)
		{
			for(j = i; j < *lenFreq; ++j)
			{
				allelePair[inc] = i; //*(allelePair + inc) = i;
				allelePair[inc + 1] = j;
				inc = inc + 2;
			}
		}
		
		double *Uset;
		int lenUset = (*x)*2;
		double *Gg;
		int lenGg;
		if(*lenT != 1)
		{
			lenGg = lenUset + *lenT;//ith possible combination known + unknown
			Gg = MALLOC(double, lenGg);
			memcpy(Gg, T, (*lenT)*sizeof(double));
			Uset = &(Gg[*lenT]);//Uset = Gg[lenT, lenT + lenUset - 1]
		}
		else//if T=0
		{
			lenGg = lenUset;
			Gg = MALLOC(double, lenUset);
			Uset = Gg;
		}
#ifdef DEBUG_TOOL
		PRINT_SET_FUNCTION("%f", Gg, lenGg);
#endif
		
		int *incrementForEachPair = MALLOC(int, *x);// vector of len x
		for(i = 0; i < *x; ++i)
		{
			incrementForEachPair[i] = 0;
		}
		
		double sumResGeno = 0.;

#ifdef DEBUG_TOOL
		int count = 0;
#endif

		do
		{
			for(i = 0; i < (*x)*2; i = i + 2)
			{
				Uset[i] = Q[allelePair[2*incrementForEachPair[i/2]]];
				Uset[i + 1] = Q[allelePair[2*incrementForEachPair[i/2] + 1]];
			}

#ifdef DEBUG_TOOL
		PRINT_SET_FUNCTION("%d", incrementForEachPair, *x);
		PRINT_SET_FUNCTION("%f", Gg, lenGg);
		PRINT_SET_FUNCTION("%f", Uset, lenUset);
#endif
		
			// for a given genotypic combination, calculate the product of the replicates probabilities.
			// genotype probability: Pelvis.
			double a = Pevid8(Uset, lenUset, Q, frequence, *lenFreq, T, *lenT, V, *lenV, *theta);

#ifdef DEBUG_TOOL
		fprintf(file, "a %f\r\n", a);
		fflush(file);
#endif

			double prodResRep = infoRepLoop2(Repliste, *lenRepliste, Gg, lenGg, prDHet, *lenHet
				, prDHom, *lenHom, *prC, Q, frequence, *lenFreq);

#ifdef DEBUG_TOOL
		fprintf(file, "prodResRep %f\r\n", prodResRep);
		fflush(file);
#endif

			sumResGeno = sumResGeno + a*prodResRep;
			
			for(i = *x - 1; i >= 0; --i)
			{
				incrementForEachPair[i]++;
				if(incrementForEachPair[i] >= nbAllelePair)
				{
#ifdef WITH_PERMUTATIONS
					incrementForEachPair[i] = 0;
#else
					if(i != 0)
					{
						for(j = i - 1; j >= 0; --j)
						{
							incrementForEachPair[i] = incrementForEachPair[j] + 1;
							if(incrementForEachPair[i] < nbAllelePair)
							{
								break;
							}
						}
					}
#endif
					continue;
				}
				else
				{
					break;
				}
			}

#ifdef DEBUG_TOOL
			count++;
#endif
		}while(i != -1);

#ifdef DEBUG_TOOL
		fprintf(file, "count %d\r\n", count);
		fflush(file);
#endif

		*sortieR = sumResGeno;
		
		free(allelePair);
		allelePair = NULL;
		free(Gg);
		Gg = NULL;
		free(incrementForEachPair);
		incrementForEachPair = NULL;	
	}
	// if x equal 0, no need to call Pevid.
	else
	{
		// if there are no unknown contributors to the sample, then just return 1 x Pr(replicate).
		double prodResRep = infoRepLoop2(Repliste, *lenRepliste, T, *lenT, prDHet, *lenHet
			, prDHom, *lenHom, *prC, Q, frequence, *lenFreq);
		*sortieR = prodResRep;
	}

	free(Q);
	Q = NULL;
	
#ifdef DEBUG_TOOL
	fprintf(file, "end debug evidenceC\n");
	fclose(file);
#endif
}
