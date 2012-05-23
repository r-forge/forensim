#include <R.h>
#include <math.h>
#include <Rmath.h>

// #define DEBUG_TOOL

#ifdef DEBUG_TOOL
	FILE *file;
#endif

double Product(double *x, int length)
{
	int i=0;
	double tmp=1.;
	while(i<length)
	{

	tmp=tmp*x[i];
	i++;
	}
	return tmp;
}

// x power y (x^y)
double Pow(double x, int y)
{
	double tmp = 1.;
	while(y)
	{
		tmp = tmp*x;
		y--;
	}
	return tmp;
}

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

void DiffIntersect(double *setInter1, int lenSetInter1
				  ,double *setInter2, int lenSetInter2
				  ,double **res, int *len)
{
	//fprintf(file, "%f %d %f %d\r\n", setInter1[0], lenSetInter1, setInter2[0], lenSetInter2);
	(*res) = (double*)malloc(lenSetInter1*sizeof(double));// assign memory to res
	*len = 0;
	int i;
	int j;
	for(i = 0; i < lenSetInter1; i++)
	{
		for(j = 0; j < lenSetInter2; j++)
		{
			if(setInter1[i] == setInter2[j])
			{
				break;
			}
		}
		if(j == lenSetInter2)
		{
			(*res)[*len]=setInter1[i];
			*len = (*len) + 1;
		}
	}
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

#define MALLOC(_type, _size) (_type*)malloc((_size)*sizeof(_type))

void infoRepC(double *R, int *lenR, double *Ggi, int *lenGgi, double *prDHet, int *lenHet, double *prDHom, int *lenHom, double *prC, char **allele, double *frequence, int *lenFreq, double *res)
{
	#ifdef DEBUG_TOOL
		file = fopen("debug2.txt","w+b");
		fprintf(file, "begin debug2\r\n");
		fflush(file);
		if(R != NULL)
		{
			fprintf(file, "res = %f, R=%f\r\n", *res, *R);
			fflush(file);
			fprintf(file, "R1 = %f, R3 =%f\r\n", R[0], R[2]);
			fflush(file);
		}
		fprintf(file, "allele = %s, allele=%s\r\n", allele[0], allele[1]);
		fflush(file);
		fprintf(file, "frequence = %f, %f\r\n", frequence[0], frequence[1]);
		fflush(file);
	#endif

	int i = 0;
	int j = 0;
	int k = 0;
	//int l = 0;
	
	double *alleleDouble = (double*)malloc((*lenFreq)*sizeof(double));
	for(i = 0; i < *lenFreq; ++i)
	{
		alleleDouble[i] = atof(allele[i]); // man atof : Conversion d'une chaîne en réel (double). 

	#ifdef DEBUG_TOOL
		fprintf(file, "alleleDouble[%d] = %f\r\n", i, alleleDouble[i]);
		fflush(file);
	#endif
	}
	
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
	
	/*double *tmpVector = MALLOC(double*, lenUnionRAndGgi);
	memset(tmpVector, 0, lenUnionRAndGgi*sizeof(double));
	int currentRow = 0;
	double proba = 0.;
	for(;;)
	{
		double product = 1.;
		for(j = 0; j < lenUnionRAndGgi; ++j)
		{
			product *= probaVector[j][tmpVector[j]];
		}
		proba += product;
		
		if(tmpVector[currentRow] < powGgi - 1)
		{
			tmpVector[currentRow]++;
		}
		else if(
		{
			currentRow++;
		}
	}
	free(tmpVector);*/
	
	*res = product;
	
#ifdef DEBUG_TOOL
    fprintf(file, "end debug2\r\n");
    fclose(file);
#endif

	//for(i = 0; i < lenUnionRAndGgi; ++i)
	//{
	//	free(probaVector[i]);
	//}
	//free(probaVector);
	free(alleleDouble);
	free(prD);
	free(UnionRAndGgi);
}


 // fonction PevidC

////// Fonction recurs2
void taballoc (double ***tab, int l1, int c1)
{
    int i;
	
#ifdef DEBUG_TOOL
	if( l1 <= 0 || c1 <= 0 )
	{
		fprintf(file, "allocation error with line %d column %d\n", l1, c1);
	}
#endif

    if ( (*tab = (double **) calloc(l1/*+1*/, sizeof(double *))) != 0) {
        for (i=0;i</*=*/l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1/*+1*/, sizeof(double))) == 0 ) {
#ifdef DEBUG_TOOL
				fprintf(file, "error no memory to alloc %d bytes\n", (c1+1)*sizeof(double));
#endif
                return;
                /*for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }*/
            }
        }
    }
#ifdef DEBUG_TOOL
	else
	{
		fprintf(file, "error no memory to alloc %d bytes\n", (l1+1)*sizeof(double *));
	}
#endif

    //**(*tab) = l1;
    //**(*tab+1) = c1;
}


//eviter le warning, free demande unchar *
void freetab (double **tab, int n)
{
    int     i;//, n;

    //n = *(*(tab));
    for (i=0;i</*=*/n;i++) {
            free((char *) *(tab+i) );
    }
    free((char *) tab);

}



void vide( int *a, int len)
{
    int i;
    for( i = 0; i < len; i++)
    {
        a[i] = 0;
    }

}

int cpt = 0;

void recurs2(int r, int c, int S, int i,int *a, double **matrice, int nbLignes)
 {
	int j, k;

        if (i == c - 1)
        {
            a[ c - 1 ] = r - S;
            for(j=0;j<c;j++)
             {
#ifdef DEBUG_TOOL
                if( cpt < 0 || cpt >= nbLignes)
                {
                    fprintf(file, "error try to reach line %d although nbligne = %d\n", cpt + 1, nbLignes);
                }
#endif
                matrice[ cpt][j]=a[j];
             }
            cpt++;
        }
        else
        {
            for (k=0; k <= r; k++)
            {
                if (S + k <= r)
                {
                    a[i] = k;

                    recurs2(r, c, S + k, i + 1, a, matrice, nbLignes);
                }
                else
                {
                    break;
                }
            }
        }


 }

void recurs_C(int r, int c, double **mat, int nbLigne, int nbCol)
{
    int *a = malloc(sizeof(int)*c);
    vide(a, c);
	
    cpt = 0;
	recurs2(r, c, 0, 0, a, mat, nbLigne);
	
    free(a);
}
 
void recurs(int *r, int *c, int *matR, int *nbLigne, int *nbCol)
{
// #ifdef DEBUG_TOOL
    // file = fopen("/home/dibule/Recurs/debug.txt","ab");
    // fprintf(file, "begin debug recurs\n");
// #endif

    double **matrice;
	taballoc(&matrice, *nbLigne, *nbCol);
	
	recurs_C(*r, *c, matrice, *nbLigne, *nbCol);

    int inc = 0;
    int l,p;
    for(p=0;p<*nbCol;p++)
	{
		for(l=0;l<*nbLigne;l++)
		{
		   //matR[j*nbLignes + k] = matrice[j][k];
			matR[inc] = matrice[l][p];
			inc=inc+1;
		}
	}

    freetab(matrice, *nbLigne);

// #ifdef DEBUG_TOOL
    // for( int i = 0; i < (*nbCol)*(*nbLigne); i++ )
    // {
        // fprintf(file, "%d ", matR[i]);
    // }
    // fprintf(file, "end debug recurs\n");
    // fclose(file);
// #endif
}







// product(x,y)


void Productby(double *x, double y, int lenx)
{
	int j;
	//avant
	// #ifdef DEBUG_TOOL
	// {
		// int i;
		// for(i = 0; i < lenx; ++i)
		// {
			// fprintf(file, "x[%d] = %e\r\n", i, x[i]);
			// fprintf(file, "vi[%d] = %e\r\n", i, vi[i]);

		// }
	// }
	// #endif
	for(j = 0; j < lenx; ++j)
	{
		x[j] = x[j]*y;
	}
	// apres
		// #ifdef DEBUG_TOOL
	// {
		// int i;
		// for(i = 0; i < lenx; ++i)
		// {
			// fprintf(file, "x[%d] = %e\r\n", i, x[i]);
			// fprintf(file, "vi[%d] = %e\r\n", i, vi[i]);

		// }
	// }
	// #endif
}


// fonction Factorial
double Fact(int n)
{
	double tmp = 1;
	while(n > 1)
	{
		tmp = tmp*n;
		--n;
	}
	return tmp;
}
// fonction Cmn: combinations with replacement
// int Cmn(int m, int n)
// {
	// return Fact(n+m-1)/(Fact(n-1)*Fact(m));
// }

//function to calcule the probability of Pr(T,U,V|H); but used as: Pr(U|T,V,H)
//-----------------------------------------------------------------------

// function to calculate the # of alleles of a given genotype  in  a stain
double* AlleleCount(double *Stain, int lenStain, double *Geno, int lenGeno )
{
// je mets les * car ce sont des vecteurs (??)
	int i=0;
	int j;
	double *res = (double*)malloc((lenStain)*sizeof(double));//create var of length lenStain
	memset(res, 0, (lenStain)*sizeof(double));

	// if(lenGeno!=0)
	// {
		// initialiser res a zero partout
		// for(k=0; k<*lenStain;k++)
		// { 
			// (*res)[k]=0;
		// }
		
		for(i=0;i<lenStain;i++)
		{
			for(j=0;j<lenGeno;j++)
			{
				if(Geno[j]==Stain[i])
				{
					res[i]=res[i]+1;

				}
			}
		
		
		}
	// }
	return res;
	
}



void PevidC(double *stain, int *lenStain, double *freq, int *lenFreq, int *x, double *T, int *lenT, double *V, int *lenV, double *theta,int *nbLigne, double *results)
{
	#ifdef DEBUG_TOOL
		file = fopen("debug.txt","wb");
		fprintf(file, "begin debug\r\n");
		fflush(file);
		

		
		fprintf(file, "stain[0]=%e, stain[1]=%e, lenStain=%d, \
		freq[0]=%e, freq[1]=%e, lenFreq=%d, x=%d, \
		T[0]=%e, T[1]=%e, lenT=%d, V[0]=%e, V[1]=%e, lenV=%d, theta=%e\r\n",
		stain[0], stain[1], *lenStain, freq[0], freq[1], *lenFreq, *x, T[0], T[1], *lenT,
		V[0], V[1], *lenV, *theta);
		fflush(file);

		// fprintf(file, "T = %e\n", *lenT);
	
	#endif
	
	
	
	double *vi= AlleleCount(stain,*lenStain,V,*lenV);
	double *ti = AlleleCount(stain,*lenStain,T,*lenT);

	
	
	// #ifdef DEBUG_TOOL
	// {
		// int i;
		// for(i = 0; i < *lenStain; ++i)
		// {
			// fprintf(file, "ti[%d] = %e\r\n", i, ti[i]);
		// }
		// for(i = 0; i < *lenStain; ++i)
		// {
			// fprintf(file, "vi[%d] = %e\r\n", i, vi[i]);
		// }
	// }
	// #endif
	
	//numbers of contributors: known, and known non-contributors
	int nT=*lenT;///2;
	int nV=*lenV;///2;
	
	// double *ptr = AlleleCount(stain,*lenStain,T,*lenT);
	
	// int condstop=;
	// double *const0;
	// ce qui suit n'est qu'un test
	// for(i=0;i<*lenStain;i++)
	// {
		// foo[i]= ti[i];
		// *(foo+i)=*(ptr+i);
	// }
	
	// calculate the constant in Curran's formula: x=nU in the original formula

	double const0=0;
	if(*x==0)
	{
		const0 = Fact(2*(*x));

	}
	else
	{
		double *tmp = (double*)malloc((2*(*x))*sizeof(double));//create var of length lenStain
		int k = 0;
		int j;
		for(j=2*nT + 2*nV;    j<=2*nT + 2*nV + 2*(*x)-1;     j++)
		{
			tmp[k]= (1 - *theta) + j * (*theta);
			k++;
		}
		const0=Fact(2 *(*x))/Product(tmp, 2*(*x));
		free(tmp);

	}
	// #ifdef DEBUG_TOOL
	// {
		
		// fprintf(file, "const0 = %e\r\n", const0);
		
	// }
	// #endif
	double ** mat = NULL;//adresse du double pointeur ca fait trois pointeurs
	//r=2x
	//c=lenStain
	// int nbLigne = Cmn(2*(*x),*lenStain);
	taballoc(&mat, *nbLigne, *lenStain);// 
		// recurs(2*(&x), &lenStain, &mat, &nbLigne, &lenStain);//modifies the matR matrix, some redundncies that will be corrected later

	recurs_C(2*(*x), *lenStain, mat, *nbLigne, *lenStain);//modifies the matR matrix, some redundncies that will be corrected later
	// void recurs(int *r, int *c, int *matR, int *nbLigne, int *nbCol)
	// r ici =2x, ce n'est pas toujours le cas
	//nblignes Cmn
	
	//double *results = (double*)malloc((nbLigne)*sizeof(double));//create var of length lenStain // prod_p <- rep(0, c)
	int l;
	int p;
	double *prod_p = (double*)malloc((*lenStain)*sizeof(double));//create var of length lenStain // prod_p <- rep(0, c)
	for(l=0;l<*nbLigne;l++)
	{	
		double productFact = 1;
		for(p=0;p<*lenStain;p++)
		{
			if(mat[l][p]<=0)
			{
				prod_p[p]= 1;
			}
			else
			{
				double *tmp2 = (double*)malloc(mat[l][p]*sizeof(double));//create var of length lenStain
				int s = 0;
				int w;
				for(w=ti[p] + vi[p];    w<=ti[p] + vi[p] + mat[l][p] -1;     w++)
				{
									
					// #ifdef DEBUG_TOOL
					// {
						
						// fprintf(file, "w=%d\r\n", w);
					// }
					// #endif
					
					
					
					tmp2[s]= (1 - *theta)*freq[p] + w * (*theta);
						s++;
				}
				// #ifdef DEBUG_TOOL
				// {
					// int i;
					// for(i = 0; i < mat[l][p]; ++i)
					// {
						// fprintf(file, "tmp2[%d] = %e\r\n", i, tmp2[i]);
					// }
				// }
				// #endif
		
				//((ti[p] + vi[i]):(ti[p] + vi[i] + ui[d, i] -1))
				prod_p[p] = Product(tmp2,mat[l][p]);
				free(tmp2);
			}
			productFact = productFact*Fact(mat[l][p]);
		}
		results[l] = Product(prod_p, *lenStain)/productFact;//prod(prod_p) / prod(factorial(ui[d, ]))
		// #ifdef DEBUG_TOOL
		// {
			
			// fprintf(file, "results = %e\r\n", results[l]);
			// int i;
			// for(i = 0; i < *lenStain; ++i)
			// {
			// fprintf(file, "prod_p[%d] = %e\r\n", i, prod_p[i]);
			// }
			
		// }
		// #endif
		
	}

	free(prod_p);
	
	Productby(results,const0,*nbLigne);

	// #ifdef DEBUG_TOOL
	// {
		// int i;
		// for(i = 0; i < *nbLigne; ++i)
		// {
			
			// fprintf(file, "results[%d] = %e\r\n", i, results[i]);
		// }
		// fprintf(file, "nbL=%d\r\n", *nbLigne);
	// }
	// #endif
	
	//free pointers
	freetab(mat, *nbLigne);
	free(ti);
	free(vi);
	//free(results);
	
	
	#ifdef DEBUG_TOOL
		fprintf(file, "z\r\n");
		fflush(file);
	#endif
	
	#ifdef DEBUG_TOOL
	// fprintf(file, "foo=%e, foo=%e, foo=%e\n",foo[0], foo[1], foo[2]);
    fprintf(file, "end debug\r\n");
    fclose(file);
	#endif
}
