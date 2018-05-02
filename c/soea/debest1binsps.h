#ifndef _DE_H_
#define _DE_H_

#define DE_F  (0.7)
#define DE_CR (0.5)

typedef struct
{
    double (*fitfun)(double*);
    double **X;
    double *fx;
    double **SP;
    double *fSP;
    uint32_t iSP;
    uint32_t *q;
    uint32_t D;
    uint32_t NP;
    uint32_t Q;
    uint32_t maxfunevals;

} debest1binsps_t;

#endif