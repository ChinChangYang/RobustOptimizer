#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "debest1binsps.h"

double sphere(double *x, uint32_t D)
{
    uint32_t i;
    double y = 3.0;

    for (i = 0; i < D; i++)
    {
        double z;

        z = x[i] - 5.0;
        y += z * z;
    }

    return y;
}

double sphere2d(double *x)
{
    return sphere(x, 2u);
}

double rastrigin(double *x, uint32_t D)
{
    uint32_t i;
    double y = 3.0;

    for (i = 0; i < D; i++)
    {
        double z;

        z = x[i] - 5.0;
        y += (z * z) - (10.0 * cos(2 * M_PI * z));
    }

    y += 10.0 * D;

    return y;
}

double rastrigin2d(double *x)
{
    return rastrigin(x, 2u);
}

double randdouble(double min, double max)
{
    double retval;
    
    retval = (double)rand();
    retval /= RAND_MAX;
    retval *= max - min;
    retval += min;

    return retval;
}

uint32_t findmin(double *fx, uint32_t NP)
{
    uint32_t i;
    uint32_t minidx;
    double minval;

    minval = fx[0];
    minidx = 0u;

    for (i = 1u; i < NP; i++)
    {
        if (fx[i] < minval)
        {
            minval = fx[i];
            minidx = i;
        }
    }

    return minidx;
}

int debest1binsps(debest1binsps_t de)
{
    uint32_t i, j, r1, r2, jrand, counteval;
    double *fu;
    double **V, **U;

    counteval = 0u;

    V = (double**)calloc(de.NP, sizeof(double*));

    for (i = 0u; i < de.NP; i++)
    {
        V[i] = (double*)calloc(de.D, sizeof(double));
    }

    U = (double**)calloc(de.NP, sizeof(double*));

    for (i = 0u; i < de.NP; i++)
    {
        U[i] = (double*)calloc(de.D, sizeof(double));
    }

    fu = (double*)calloc(de.NP, sizeof(double));

    while (1)
    {
        uint32_t ibestX, ibestSP;

        /* Termination conditions */
        if (counteval > (de.maxfunevals - de.NP))
        {
            break;
        }

        /* Mutation */
        ibestX = findmin(de.fx, de.NP);
        ibestSP = findmin(de.fSP, de.NP);

        for (i = 0u; i < de.NP; i++)
        {
            /* Generate r1 */
            r1 = (uint32_t)(randdouble(0.0, de.NP - 1.0));

            while (i == r1)
            {
                r1 = (uint32_t)(randdouble(0.0, de.NP - 1.0));
            }

            /* Generate r2 */
            r2 = (uint32_t)(randdouble(0.0, de.NP - 1.0));

            while ((i == r2) || (r1 == r2))
            {
                r2 = (uint32_t)(randdouble(0.0, de.NP - 1.0));
            }

            if (de.q[i] <= de.Q)
            {
                for (j = 0u; j < de.D; j++)
                {
                    V[i][j] = de.X[ibestX][j] + DE_F * (de.X[r1][j] - de.X[r2][j]);
                }
            }
            else
            {
                for (j = 0u; j < de.D; j++)
                {
                    V[i][j] = de.SP[ibestSP][j] + DE_F * (de.SP[r1][j] - de.SP[r2][j]);
                }
            }
        } /* Mutation */

        for (i = 0u; i < de.NP; i++)
        {
            /* Binomial Crossover */
            jrand = (uint32_t)(randdouble(0.0, de.D - 1.0));

            if (de.q[i] <= de.Q)
            {
                for (j = 0u; j < de.D; j++)
                {
                    if ((randdouble(0.0, 1.0) < DE_CR)
                        || (j == jrand))
                    {
                        U[i][j] = V[i][j];
                    }
                    else
                    {
                        U[i][j] = de.X[i][j];
                    }
                }
            }
            else
            {
                for (j = 0u; j < de.D; j++)
                {
                    if ((randdouble(0.0, 1.0) < DE_CR)
                        || (j == jrand))
                    {
                        U[i][j] = V[i][j];
                    }
                    else
                    {
                        U[i][j] = de.SP[i][j];
                    }
                }
            }
        } /* Crossover */

        /* Evaluation */
        for (i = 0u; i < de.NP; i++)
        {
            fu[i] = de.fitfun(U[i]);
            counteval++;
        }

        /* Selection */
        for (i = 0u; i < de.NP; i++)
        {
            if (fu[i] < de.fx[i])
            {
                for (j = 0u; j < de.D; j++)
                {
                    de.X[i][j] = U[i][j];
                    de.SP[de.iSP][j] = U[i][j];
                }

                de.fx[i] = fu[i];
                de.fSP[de.iSP] = fu[i];
                de.iSP = de.iSP % de.NP;
                de.q[i] = 0u;
            }
            else
            {
                de.q[i]++;
            }
        }

    } /* while (true) */

    for (i = 0u; i < de.D; i++)
    {
        free(V[i]);
        free(U[i]);
    }

    free(V);
    free(U);
    free(fu);

    return 0;
}

int main()
{
    uint32_t i, j;
    debest1binsps_t de;

    de.D = 2u;
    de.NP = 10u;
    de.maxfunevals = 1000u;
    de.Q = 64u;
    de.iSP = 0u;
    de.fitfun = &rastrigin2d;
    de.fx = (double*)calloc(de.NP, sizeof(double));
    de.fSP = (double*)calloc(de.NP, sizeof(double));
    de.q = (uint32_t*)calloc(de.NP, sizeof(uint32_t));
    
    de.X = (double**)calloc(de.NP, sizeof(double*));
    de.SP = (double**)calloc(de.NP, sizeof(double*));

    for (i = 0u; i < de.NP; i++)
    {
        de.X[i] = (double*)calloc(de.D, sizeof(double));
        de.SP[i] = (double*)calloc(de.D, sizeof(double));

        for (j = 0u; j < de.D; j++)
        {
            de.X[i][j] = randdouble(-100.0, 100.0);
            de.SP[i][j] = de.X[i][j];
        }

        de.fx[i] = de.fitfun(de.X[i]);
        de.fSP[i] = de.fitfun(de.SP[i]);
    }

    debest1binsps(de);

    for (i = 0u; i < de.NP; i++)
    {
        printf("fx[%u]=%.4E\n", i, de.fx[i]);
    }

    for (i = 0u; i < de.D; i++)
    {
        free(de.X[i]);
        free(de.SP[i]);
    }

    free(de.X);
    free(de.SP);
    free(de.fx);
    free(de.fSP);
    free(de.q);

    return 0;
}
