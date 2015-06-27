/*
 * CEC15 Test Function Suite
 * Jane Jing Liang (email: liangjing@zzu.edu.cn)
 * Nov. 16th 2014
 * Reference: J. J. Liang, B. Y. Qu, P. N. Suganthan, Q. Chen, "Problem Definitions and Evaluation Criteria for the CEC 2015 Competition
 * on Learning-based Real-Parameter Single Objective Optimization",Technical Report201411A, Computational Intelligence
 * Laboratory, Zhengzhou University, Zhengzhou, China and Technical Report, Nanyang Technological University, Singapore, 2014.
 * 1. Run "mex -setup" and choose a proper complier (Microsoft Visual C++ 6.0 is prefered)
 * Would you like mex to locate installed compilers [y]/n? y
 * 2. Run the following command in Matlab window:
 * mex cec15_func.cpp -DWINDOWS
 * 3. Then you can use the test functions as the following example:
 * f = cec15_func(x,func_num);
 * Here x is a D*pop_size matrix.
 */
//#include <WINDOWS.H>
#include <iostream>
#include <fstream>
#include <math.h>
#include <mex.h>

double *OShift=NULL;
double *M=NULL;
double *y=NULL;
double *z=NULL;
double *x_bound=NULL;
double *bias =NULL;
int ini_flag=0,n_flag=-1,func_flag;
int*SS=NULL;


#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void exitFcn()
{
    if (M!=NULL)
        mxFree(M);
    if (bias!=NULL)
        mxFree(bias);
    if (SS !=NULL)
        mxFree(SS);
    if (OShift!=NULL)
        mxFree(OShift);
    if (y!=NULL)
        mxFree(y);
    if (z!=NULL)
        mxFree(z);
    if (x_bound!=NULL)
        mxFree(x_bound);
}
void makeMemory()
{
    
    M=(double*)mxCalloc(10*100*100,sizeof(double));
    bias=(double*)mxCalloc(10,sizeof(double));
    SS=(int *)mxCalloc(100*10,sizeof(int));
    OShift=(double *)mxCalloc(100*10,sizeof(double));
    y=(double *)mxCalloc(100, sizeof(double));
    z=(double *)mxCalloc(100, sizeof(double));
    x_bound=(double *)mxCalloc(100, sizeof(double));
    for (int i=0; i<100; i++)
        x_bound[i]=100.0;
    mexMakeMemoryPersistent(M);
    mexMakeMemoryPersistent(bias);
    mexMakeMemoryPersistent(SS);
    mexMakeMemoryPersistent(OShift);
    mexMakeMemoryPersistent(y);
    mexMakeMemoryPersistent(z);
    mexMakeMemoryPersistent(x_bound);
    
}



void sphere_func (double *, double *, int , double *,double *, int, int); /* Sphere */
void ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int , double *,double *, int, int); /* Bent_Cigar */
void discus_func(double *, double *, int , double *,double *, int, int);  /* Discus */
void dif_powers_func(double *, double *, int , double *,double *, int, int);  /* Different Powers */
void rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void schaffer_F7_func (double *, double *, int , double *,double *, int, int); /* Schwefel's F7 */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void weierstrass_func (double *, double *, int , double *,double *, int, int); /* Weierstrass's  */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void katsuura_func (double *, double *, int , double *,double *, int, int); /* Katsuura */
void bi_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Lunacek Bi_rastrigin */
void grie_rosen_func (double *, double *, int , double *,double *, int, int); /* Griewank-Rosenbrock  */
void escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer¡¯s F6  */
void step_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Noncontinuous Rastrigin's  */
void happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void hgbat_func (double *, double *, int , double *,double *, int, int); /* HGBat  */

void hf01 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 1 */
void hf02 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 2 */
void hf03 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 3 */

void cf01 (double *, double *, int , double *,double *,double *, int); /* Composition Function 1 */
void cf02 (double *, double *, int , double *,double *, int *,double *, int); /* Composition Function 2 */
void cf03 (double *, double *, int , double *,double *, double *,int); /* Composition Function 3 */
void cf04 (double *, double *, int , double *,double *, double *,int); /* Composition Function 4 */
void cf05 (double *, double *, int , double *,double *, int *,double *, int); /* Composition Function 5 */
void cf06 (double *, double *, int , double *,double *, double *,int); /* Composition Function 6 */
void cf07 (double *, double *, int , double *,double *, double *,int); /* Composition Function 7 */

void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void cf_cal(double *, double *, int, double *,double *,double *,double *,int);
void cec15_test_func(double *, double *,int,int,int);

#include <exception>
#include <limits>

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    
    int  m, n,func_num;
    double  *f, *x;
    
    mexAtExit(&exitFcn);
    if (nrhs < 2)
    {
        mexPrintf ("usage: f = cec15_func(x, func_num);\n");
        mexErrMsgTxt ("example: f= cec15_func([3.3253000e+000, -1.2835000e+000]', 1);");
    }
    n = (int)mxGetM (prhs[0]);
    m = (int)mxGetN (prhs[0]);
    x = (double*)mxGetPr (prhs[0]);
    func_num= (int)*mxGetPr (prhs[1]);
    if (func_num>15)
    {
        mexPrintf ("usage: f = cec15_func(x, func_num);\n");
        mexErrMsgTxt ("Error: There are only 15 test functions in this test suite!");
    }
    
    plhs[0] = mxCreateDoubleMatrix (1, m, mxREAL);
    f = (double*)mxGetPr (plhs[0]);
    cec15_test_func(&x[0], &f[0], n,m,func_num);
    
}


void cec15_test_func(double *x, double *f, int nx, int mx,int func_num)
{
    if (M==NULL)
        makeMemory();
    int cf_num=10,i,j;
    std::ifstream fpt;
    
    if (ini_flag==1)
    {
        if ((n_flag!=nx)||(func_flag!=func_num))
        {
            ini_flag=0;
        }
    }
    
    if (ini_flag==0)
    {
        char FileName[256];
        if (!(nx==2||nx==10||nx==30||nx==50||nx==100))
        {
            mexErrMsgTxt("\nError: Test functions are only defined for D=2,10,30,50,100.\n");
        }
        if (nx==2&&((func_num>=6&&func_num<=8)||(func_num==10)||(func_num==13)))
        {
            mexErrMsgTxt("\nError: hf01,hf02,hf03,cf02&cf05 are NOT defined for D=2.\n");
        }
        
        /* Load Matrix M*/
        sprintf(FileName, "function/soea/cec2015/input_data/M_%d_D%d.txt", func_num,nx);
        fpt.open(FileName);
        if (!fpt.good())
        {
            mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
        }
        if (func_num<9)
        {
            if (M==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for (i=0; i<nx*nx; i++)
            {
                fpt >> M[i];
            }
        }
        else
        {
            if (M==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for (i=0; i<cf_num*nx*nx; i++)
            {
                fpt >> M[i];
            }
        }
        fpt.close();
        
        /* Load Bias_value bias*/
        if (func_num>=9)
        {
            sprintf(FileName, "function/soea/cec2015/input_data/bias_%d.txt", func_num);
            fpt.open(FileName);
            if (!fpt.good())
            {
                mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
            }
            
            if (bias==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for (i=0; i<cf_num; i++)
            {
                fpt >> bias[i];
            }
            fpt.close();
        }
        
        /* Load shift_data */
        sprintf(FileName, "function/soea/cec2015/input_data/shift_data_%d.txt", func_num);
        fpt.open(FileName);
        if (!fpt.good())
        {
            mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
        }
        
        if (func_num<9)
        {
            if (OShift==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for(i=0;i<nx;i++)
            {
                fpt >> OShift[i];
            }
        }
        else
        {
            if (OShift==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for(i=0;i<cf_num*nx;i++)
            {
                fpt >> OShift[i];
                // ignore the rest of the line to the next line of data
                if (((i+1) % nx) == 0) {
                    fpt.ignore( std::numeric_limits<std::streamsize>::max(), '\n');
                }
            }
            
        }
        fpt.close();
        
        
        /* Load Shuffle_data */
        
        if (func_num>=6&&func_num<=8)
        {
            sprintf(FileName, "function/soea/cec2015/input_data/shuffle_data_%d_D%d.txt", func_num, nx);
            fpt.open(FileName);
            if (!fpt.good())
            {
                mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
            }
            if (SS==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for(i=0;i<nx;i++)
            {
                fpt >> SS[i];
            }
            fpt.close();
        }
        else if (func_num==10||func_num==13)
        {
            sprintf(FileName, "function/soea/cec2015/input_data/shuffle_data_%d_D%d.txt", func_num, nx);
            fpt.open(FileName);
            if (!fpt.good())
            {
                mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
            }
            if (SS==NULL)
                mexErrMsgTxt("\nError: there is insufficient memory available!\n");
            for(i=0;i<nx*cf_num;i++)
            {
                fpt >> SS[i];
            }
            fpt.close();
        }
        
        
        n_flag=nx;
        func_flag=func_num;
        ini_flag=1;
    }
    
    
    for (i = 0; i < mx; i++)
    {
        switch(func_num)
        {
            case 1:
                ellips_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
                f[i]+=100.0;
                break;
            case 2:
                bent_cigar_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
                f[i]+=200.0;
                break;
            case 3:
                ackley_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
                f[i]+=300.0;
                break;
            case 4:
                rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
                f[i]+=400.0;
                break;
            case 5:
                schwefel_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
                f[i]+=500.0;
                break;
            case 6:
                hf01(&x[i*nx],&f[i],nx,OShift,M,SS,1,1);
                f[i]+=600.0;
                break;
            case 7:
                hf02(&x[i*nx],&f[i],nx,OShift,M,SS,1,1);
                f[i]+=700.0;
                break;
            case 8:
                hf03(&x[i*nx],&f[i],nx,OShift,M,SS,1,1);
                f[i]+=800.0;
                break;
            case 9:
                cf01(&x[i*nx],&f[i],nx,OShift,M,bias,1);
                f[i]+=900.0;
                break;
            case 10:
                cf02(&x[i*nx],&f[i],nx,OShift,M,SS,bias,1);
                f[i]+=1000.0;
                break;
            case 11:
                cf03(&x[i*nx],&f[i],nx,OShift,M,bias,1);
                f[i]+=1100.0;
                break;
            case 12:
                cf04(&x[i*nx],&f[i],nx,OShift,M,bias,1);
                f[i]+=1200.0;
                break;
            case 13:
                cf05(&x[i*nx],&f[i],nx,OShift,M,SS,bias,1);
                f[i]+=1300.0;
                break;
            case 14:
                cf06(&x[i*nx],&f[i],nx,OShift,M,bias,1);
                f[i]+=1400.0;
                break;
            case 15:
                cf07(&x[i*nx],&f[i],nx,OShift,M,bias,1);
                f[i]+=1500.0;
                break;
                
            default:
                printf("\nError: There are only 30 test functions in this test suite!\n");
                f[i] = 0.0;
                break;
        }
        
    }
    
}
void sphere_func (double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag) /* Sphere */
{
    int i;
    f[0] = 0.0;
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    for (i=0; i<nx; i++)
    {
        f[0] += z[i]*z[i];
    }
    
}



void ellips_func (double *x, double *f, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Ellipsoidal */
{
    int i;
    f[0] = 0.0;
    sr_func (x, z, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
    for (i=0; i<nx; i++)
    {
        f[0] += pow(10.0,6.0*i/(nx-1))*z[i]*z[i];
    }
}

void bent_cigar_func (double *x, double *f, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Bent_Cigar */
{
    int i;
    sr_func (x, z, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
    f[0] = z[0]*z[0];
    for (i=1; i<nx; i++)
    {
        f[0] += pow(10.0,6.0)*z[i]*z[i];
    }
    
}

void discus_func (double *x, double *f, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Discus */
{
    int i;
    sr_func (x, z, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
    f[0] = pow(10.0,6.0)*z[0]*z[0];
    for (i=1; i<nx; i++)
    {
        f[0] += z[i]*z[i];
    }
}

void dif_powers_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Different Powers */
{
    int i;
    f[0] = 0.0;
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        f[0] += pow(fabs(z[i]),2+4*i/(nx-1));
    }
    f[0]=pow(f[0],0.5);
}


void rosenbrock_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rosenbrock's */
{
    int i;
    double tmp1,tmp2;
    f[0] = 0.0;
    sr_func (x, z, nx, Os, Mr, 2.048/100.0, s_flag, r_flag); /* shift and rotate */
    z[0] += 1.0;//shift to orgin
    for (i=0; i<nx-1; i++)
    {
        z[i+1] += 1.0;//shift to orgin
        tmp1=z[i]*z[i]-z[i+1];
        tmp2=z[i]-1.0;
        f[0] += 100.0*tmp1*tmp1 +tmp2*tmp2;
    }
}

void schaffer_F7_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's 1.2  */
{
    int i;
    double tmp;
    f[0] = 0.0;
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    for (i=0; i<nx-1; i++)
    {
        z[i]=pow(y[i]*y[i]+y[i+1]*y[i+1],0.5);
        tmp=sin(50.0*pow(z[i],0.2));
        f[0] += pow(z[i],0.5)+pow(z[i],0.5)*tmp*tmp ;
    }
    f[0] = f[0]*f[0]/(nx-1)/(nx-1);
}

void ackley_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;
    
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        sum1 += z[i]*z[i];
        sum2 += cos(2.0*PI*z[i]);
    }
    sum1 = -0.2*sqrt(sum1/nx);
    sum2 /= nx;
    f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}


void weierstrass_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    f[0] = 0.0;
    
    sr_func (x, z, nx, Os, Mr, 0.5/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        sum = 0.0;
        sum2 = 0.0;
        for (j=0; j<=k_max; j++)
        {
            sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(z[i]+0.5));
            sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
        }
        f[0] += sum;
    }
    f[0] -= nx*sum2;
}


void griewank_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;
    
    sr_func (x, z, nx, Os, Mr, 600.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        s += z[i]*z[i];
        p *= cos(z[i]/sqrt(1.0+i));
    }
    f[0] = 1.0 + s/4000.0 - p;
}

void rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rastrigin's  */
{
    int i;
    f[0] = 0.0;
    
    sr_func (x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        f[0] += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
    }
}

void step_rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
{
    int i;
    f[0]=0.0;
    for (i=0; i<nx; i++)
    {
        if (fabs(y[i]-Os[i])>0.5)
            y[i]=Os[i]+floor(2*(y[i]-Os[i])+0.5)/2;
    }
    
    sr_func (x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        f[0] += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
    }
}

void schwefel_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's  */
{
    int i;
    double tmp;
    f[0]=0.0;
    
    sr_func (x, z, nx, Os, Mr, 1000.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        z[i] += 4.209687462275036e+002;
        if (z[i]>500)
        {
            f[0]-=(500.0-fmod(z[i],500))*sin(pow(500.0-fmod(z[i],500),0.5));
            tmp=(z[i]-500.0)/100;
            f[0]+= tmp*tmp/nx;
        }
        else if (z[i]<-500)
        {
            f[0]-=(-500.0+fmod(fabs(z[i]),500))*sin(pow(500.0-fmod(fabs(z[i]),500),0.5));
            tmp=(z[i]+500.0)/100;
            f[0]+= tmp*tmp/nx;
        }
        else
            f[0]-=z[i]*sin(pow(fabs(z[i]),0.5));
    }
    f[0] +=4.189828872724338e+002*nx;
    
}

void katsuura_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Katsuura  */
{
    int i,j;
    double temp,tmp1,tmp2,tmp3;
    f[0]=1.0;
    tmp3=pow(1.0*nx,1.2);
    
    sr_func (x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        temp=0.0;
        for (j=1; j<=32; j++)
        {
            tmp1=pow(2.0,j);
            tmp2=tmp1*z[i];
            temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
        }
        f[0] *= pow(1.0+(i+1)*temp,10.0/tmp3);
    }
    tmp1=10.0/nx/nx;
    f[0]=f[0]*tmp1-tmp1;
    
}

void bi_rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Lunacek Bi_rastrigin Function */
{
    int i;
    double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
    double *tmpx;
    tmpx=(double *)malloc(sizeof(double)  *  nx);
    s=1.0-1.0/(2.0*pow(nx+20.0,0.5)-8.2);
    mu1=-pow((mu0*mu0-d)/s,0.5);
    
    if (s_flag==1)
        shiftfunc(x, y, nx, Os);
    else
    {
        for (i=0; i<nx; i++)//shrink to the orginal search range
        {
            y[i] = x[i];
        }
    }
    for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i] *= 10.0/100.0;
    }
    
    for (i = 0; i < nx; i++)
    {
        tmpx[i]=2*y[i];
        if (Os[i] < 0.0)
            tmpx[i] *= -1.;
    }
    for (i=0; i<nx; i++)
    {
        z[i]=tmpx[i];
        tmpx[i] += mu0;
    }
    tmp1=0.0;tmp2=0.0;
    for (i=0; i<nx; i++)
    {
        tmp = tmpx[i]-mu0;
        tmp1 += tmp*tmp;
        tmp = tmpx[i]-mu1;
        tmp2 += tmp*tmp;
    }
    tmp2 *= s;
    tmp2 += d*nx;
    tmp=0.0;
    
    if (r_flag==1)
    {
        rotatefunc(z, y, nx, Mr);
        for (i=0; i<nx; i++)
        {
            tmp+=cos(2.0*PI*y[i]);
        }
        if(tmp1<tmp2)
            f[0] = tmp1;
        else
            f[0] = tmp2;
        f[0] += 10.0*(nx-tmp);
    }
    else
    {
        for (i=0; i<nx; i++)
        {
            tmp+=cos(2.0*PI*z[i]);
        }
        if(tmp1<tmp2)
            f[0] = tmp1;
        else
            f[0] = tmp2;
        f[0] += 10.0*(nx-tmp);
    }
    
    free(tmpx);
}

void grie_rosen_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank-Rosenbrock  */
{
    int i;
    double temp,tmp1,tmp2;
    f[0]=0.0;
    
    sr_func (x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    z[0] += 1.0;//shift to orgin
    for (i=0; i<nx-1; i++)
    {
        z[i+1] += 1.0;//shift to orgin
        tmp1 = z[i]*z[i]-z[i+1];
        tmp2 = z[i]-1.0;
        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
        f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    tmp1 = z[nx-1]*z[nx-1]-z[0];
    tmp2 = z[nx-1]-1.0;
    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
    f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0 ;
}


void escaffer6_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Expanded Scaffer¡¯s F6  */
{
    int i;
    double temp1, temp2;
    
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    f[0] = 0.0;
    for (i=0; i<nx-1; i++)
    {
        temp1 = sin(sqrt(z[i]*z[i]+z[i+1]*z[i+1]));
        temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z[i]*z[i]+z[i+1]*z[i+1]);
        f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z[nx-1]*z[nx-1]+z[0]*z[0]));
    temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z[nx-1]*z[nx-1]+z[0]*z[0]);
    f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
}

void happycat_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha,r2,sum_z;
    alpha=1.0/8.0;
    
    sr_func (x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    r2 = 0.0;
    sum_z=0.0;
    for (i=0; i<nx; i++)
    {
        z[i]=z[i]-1.0;//shift to orgin
        r2 += z[i]*z[i];
        sum_z += z[i];
    }
    f[0]=pow(fabs(r2-nx),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
}

void hgbat_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha,r2,sum_z;
    alpha=1.0/4.0;
    
    sr_func (x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */
    
    r2 = 0.0;
    sum_z=0.0;
    for (i=0; i<nx; i++)
    {
        z[i]=z[i]-1.0;//shift to orgin
        r2 += z[i]*z[i];
        sum_z += z[i];
    }
    f[0]=pow(fabs(pow(r2,2.0)-pow(sum_z,2.0)),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
}

void hf01 (double *x, double *f, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 1 */
{
    int i,tmp,cf_num=3;
    double fit[3];
    int G[3],G_nx[3];
    double Gp[3]={0.3,0.3,0.4};
    
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_nx[i] = ceil(Gp[i]*nx);
        tmp += G_nx[i];
    }
    G_nx[cf_num-1]=nx-tmp;
    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_nx[i-1];
    }
    
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    schwefel_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=1;
    rastrigin_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=2;
    ellips_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    f[0]=0.0;
    for(i=0;i<cf_num;i++)
    {
        f[0] += fit[i];
    }
}


void hf02 (double *x, double *f, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 2 */
{
    int i,tmp,cf_num=4;
    double fit[4];
    int G[4],G_nx[4];
    double Gp[4]={0.2,0.2,0.3,0.3};
    
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_nx[i] = ceil(Gp[i]*nx);
        tmp += G_nx[i];
    }
    G_nx[cf_num-1]=nx-tmp;
    
    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_nx[i-1];
    }
    
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    griewank_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=1;
    weierstrass_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=2;
    rosenbrock_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=3;
    escaffer6_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    
    f[0]=0.0;
    for(i=0;i<cf_num;i++)
    {
        f[0] += fit[i];
    }
}


void hf03 (double *x, double *f, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 3 */
{
    int i,tmp,cf_num=5;
    double fit[5];
    int G[5],G_nx[5];
    double Gp[5]={0.1,0.2,0.2,0.2,0.3};
    
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_nx[i] = ceil(Gp[i]*nx);
        tmp += G_nx[i];
    }
    G_nx[cf_num-1]=nx-tmp;
    
    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_nx[i-1];
    }
    
    sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
    
    for (i=0; i<nx; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    escaffer6_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=1;
    hgbat_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=2;
    rosenbrock_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=3;
    schwefel_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    i=4;
    ellips_func(&y[G[i]],&fit[i],G_nx[i],Os,Mr,0,0);
    
    f[0]=0.0;
    for(i=0;i<cf_num;i++)
    {
        f[0] += fit[i];
    }
}

void cf01 (double *x, double *f, int nx, double *Os,double *Mr,double *bias, int r_flag) /* Composition Function 1 */
{
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {20,20,20};
    
    i=0;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,0);
    i=1;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    i=2;
    hgbat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf02 (double *x, double *f, int nx, double *Os,double *Mr,int *SS,double *bias,int r_flag) /* Composition Function 2 */
{
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {10,30,50};
    i=0;
    hf01(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],&SS[i*nx],1,r_flag);
    i=1;
    hf02(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],&SS[i*nx],1,r_flag);
    i=2;
    hf03(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],&SS[i*nx],1,r_flag);
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf03 (double *x, double *f, int nx, double *Os,double *Mr,double *bias,int r_flag) /* Composition Function 3 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,10,10,20,20};
    i=0;
    hgbat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1000.0;
    i=1;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    i=2;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/4.0e+3;
    i=3;
    weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/400.0;
    i=4;
    ellips_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+10;
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf04 (double *x, double *f, int nx, double *Os,double *Mr,double *bias,int r_flag) /* Composition Function 4 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,20,20,30,30};
    i=0;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/4.0e+3;
    i=1;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    i=2;
    ellips_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+10;
    i=3;
    escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=fit[i]*10.0;
    i=4;
    happycat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf05 (double *x, double *f, int nx, double *Os,double *Mr,int *SS,double *bias,int r_flag) /* Composition Function 5 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,10,10,20,20};
    i=0;
    hf03(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],&SS[i*nx],1,r_flag);
    i=1;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    i=2;
    hf01(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],&SS[i*nx],1,r_flag);
    i=3;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/4.0e+3;
    i=4;
    escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=fit[i]*10.0;
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf06 (double *x, double *f, int nx, double *Os,double *Mr,double *bias,int r_flag) /* Composition Function 6 */
{
    int i,cf_num=7;
    double fit[7];
    double delta[7] = {10,20,30,40,50,50,50};
    i=0;
    happycat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    i=1;
    grie_rosen_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/4.0e+3;
    i=2;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/4.0e+3;
    i=3;
    escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=fit[i]*10.0;
    i=4;
    ellips_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+10;
    i=5;
    bent_cigar_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+10;
    i=6;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=10000.0*fit[i]/1.0e+3;
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf07 (double *x, double *f, int nx, double *Os,double *Mr,double *bias,int r_flag) /* Composition Function 7 */
{
    int i,cf_num=10;
    double fit[10];
    double delta[10] = {10,10,20,20,30,30,40, 40, 50, 50};
    i=0;
    rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1.0e+3;
    i=1;
    weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/400.0;
    i=2;
    happycat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1.0e+3;
    i=3;
    schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/4.0e+3;
    i=4;
    rosenbrock_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1.0e+5;
    i=5;
    hgbat_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1000.0;
    i=6;
    katsuura_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1.0e+7;
    i=7;
    escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=fit[i]*10.0;
    i=8;
    grie_rosen_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/4.0e+3;
    i=9;
    ackley_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
    fit[i]=100.0*fit[i]/1.0e+5;
    cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}


void shiftfunc (double *x, double *xshift, int nx,double *Os)
{
    int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc (double *x, double *xrot, int nx,double *Mr)
{
    int i,j;
    for (i=0; i<nx; i++)
    {
        xrot[i]=0;
        for (j=0; j<nx; j++)
        {
            xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
        }
    }
}

void sr_func (double *x, double *sr_x, int nx, double *Os,double *Mr, double sh_rate, int s_flag,int r_flag) /* shift and rotate */
{
    int i;
    if (s_flag==1)
    {
        if (r_flag==1)
        {
            shiftfunc(x, y, nx, Os);
            for (i=0; i<nx; i++)//shrink to the orginal search range
            {
                y[i]=y[i]*sh_rate;
            }
            rotatefunc(y, sr_x, nx, Mr);
        }
        else
        {
            shiftfunc(x, sr_x, nx, Os);
            for (i=0; i<nx; i++)//shrink to the orginal search range
            {
                sr_x[i]=sr_x[i]*sh_rate;
            }
        }
    }
    else
    {
        
        if (r_flag==1)
        {
            for (i=0; i<nx; i++)//shrink to the orginal search range
            {
                y[i]=x[i]*sh_rate;
            }
            rotatefunc(y, sr_x, nx, Mr);
        }
        else
            for (i=0; i<nx; i++)
            {
                for (i=0; i<nx; i++)//shrink to the orginal search range
                {
                    sr_x[i]=x[i]*sh_rate;
                }
            }
    }
}

void cf_cal(double *x, double *f, int nx, double *Os,double * delta,double * bias,double * fit, int cf_num)
{
    int i,j;
    //double *w;
    double w_max=0,w_sum=0;
    //w=(double *)malloc(cf_num * sizeof(double));
    double w[10];
    for (i=0; i<cf_num; i++)
    {
        fit[i]+=bias[i];
        w[i]=0;
        for (j=0; j<nx; j++)
        {
            w[i]+=pow(x[j]-Os[i*nx+j],2.0);
        }
        if (w[i]!=0)
            w[i]=pow(1.0/w[i],0.5)*exp(-w[i]/2.0/nx/pow(delta[i],2.0));
        else
            w[i]=INF;
        if (w[i]>w_max)
            w_max=w[i];
    }
    
    for (i=0; i<cf_num; i++)
    {
        w_sum=w_sum+w[i];
    }
    if(w_max==0)
    {
        for (i=0; i<cf_num; i++)
            w[i]=1;
        w_sum=cf_num;
    }
    f[0] = 0.0;
    for (i=0; i<cf_num; i++)
    {
        f[0]=f[0]+w[i]/w_sum*fit[i];
    }
    //free(w);
}

