#include <R.h>
#include <Rmath.h>


void taballoc (double ***tab, int l1, int c1)
{
    int i, j;

    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
        for (i=0;i<=l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
                return;
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
            }
        }
    }

    **(*tab) = l1;
    **(*tab+1) = c1;

}



void freetab (double **tab)
{
    int     i, n;

    n = *(*(tab));
    for (i=0;i<=n;i++) {
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

void recurs(int *r, int *c, int *matR, int *nbLigne, int *nbCol)
{
    double **matrice;
	taballoc(&matrice, *nbLigne, *nbCol);
    int *a = malloc(sizeof(int)*(*c));
    vide(a, *c);

    cpt = 0;
	recurs2(*r, *c, 0, 0, a, matrice, *nbLigne);

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

    freetab(matrice);
    free(a);
}








