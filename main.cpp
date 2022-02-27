#include <iostream>
#include <math.h>
#include "libstats.h"

using namespace std;

static double k = 1.0;

void Therm_Runge_Kutta(double *x,double *y,double *z,double h);
void Runge_Kutta(double *x,double *y,double *z,double h,int temporalDelay);
double RosslerX(double y,double z);
double RosslerY_Delay(double x,double y,double yDelay);
double RosslerY(double x,double y);
double RosslerZ(double x,double z);
double Distance(double x1,double x2,double y1,double y2,double z1,double z2);

int main(int argc, const char * argv[]) {
    unsigned long t0,t1;
    t0 = clock();
    
    int i,j,m,NIntegration,measuInterval,temporalDelay;
    double x1,z1,yaux1,h,t,F;
    temporalDelay = 17500;
    double *y1;
    y1 = (double*)malloc(temporalDelay*sizeof(double));

//    //Computation of the Lyanpunov exponents
//    double Lyapunov[1000],av,error;
//    double d0,d;
//    double x2,z2;
//    double *y2;
//    y2 = (double*)malloc(temporalDelay*sizeof(double));
    
    h = 0.001;
    t = 0;
    NIntegration = 100000;
    measuInterval = 20;
    
    FILE *fout;
    fout = fopen("/Users/danixala5/Documents/Master/Sistemas dinámicos y caos/Control_Caos/Positions.txt", "wt");

//    for(m=0;m<100;m++)
    {
        yaux1 = -2.147;
        x1 =  2.078;
        z1 = 27;
        for(i=0;i<temporalDelay;i++){
            Therm_Runge_Kutta(&x1, &yaux1, &z1, h);
            t+= h;
            y1[i] = yaux1;
//            y2[i] = yaux1;
        }
//        x2 = x1+pow(10,-8);
//        y2[temporalDelay-1] = y1[temporalDelay-1]+pow(10,-8);
//        z2 = z1+pow(10,-8);
//        d0 = Distance(x1,x2,y1[temporalDelay-1],y2[temporalDelay-1],z1,z2);
//        d = Distance(x1,x2,y1[temporalDelay-1],y2[temporalDelay-1],z1,z2);
//        for(j=0;j<1000;j++)
        {
//            y2[temporalDelay-1] = y1[temporalDelay-1]+(y2[temporalDelay-1]-y1[temporalDelay-1])*(d0/d);
//            x2 = x1+(x2-x1)*(d0/d);
//            z2 = z1+(z2-z1)*(d0/d);
            
            for(i=0;i<NIntegration;i++){
                Runge_Kutta(&x1, y1, &z1, h, temporalDelay);
//                Runge_Kutta(&x2, y2, &z2, h, temporalDelay);
                t+= h;
                if(i%measuInterval==0){
                    F = k*(y1[0]-y1[temporalDelay-1]);
                    if(F>0.2)
                        F = 0.2;
                    if(F<-0.2)
                        F = -0.2;
                    fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\n",t,x1,y1[temporalDelay-1],z1,F);
                }
            }
//            d = Distance(x1,x2,y1[temporalDelay-1],y2[temporalDelay-1],z1,z2);
//            Lyapunov[j] = (1.0/(NIntegration*h))*log(abs(d/d0));
        }
//        av = AverageDouble(Lyapunov, 1000);
//        error = DeviationDouble(VarianceDouble(Lyapunov, 1000, av),1000);
//        fprintf(fout,"%lf\t%lf\t%lf\n",k,av,error);
//        k+= 0.03;
//        cout << k << endl;
    }
    
    fclose(fout);
    
    t1 = clock();
    cout << "tiempo de ejecución=" << (double(t1-t0)/CLOCKS_PER_SEC) << "\n";
    
    return 0;
}

void Therm_Runge_Kutta(double *x,double *y,double *z,double h)
{
    double ky1,kx1,kz1,kx2,ky2,kz2,kx3,ky3,kz3,kx4,ky4,kz4;
    
    kx1 = RosslerX(*y, *z);
    ky1 = RosslerY(*x, *y);
    kz1 = RosslerZ(*x, *z);
    
    kx2 = RosslerX(*y+0.5*h, *z+0.5*h);
    ky2 = RosslerY(*x+0.5*h, *y+0.5*h*ky1);
    kz2 = RosslerZ(*x+0.5*h, *z+0.5*h*kz1);
    
    kx3 = RosslerX(*y+0.5*h, *z+0.5*h);
    ky3 = RosslerY(*x+0.5*h, *y+0.5*h*ky2);
    kz3 = RosslerZ(*x+0.5*h, *z+0.5*h*kz2);
    
    kx4 = RosslerX(*y+h, *z+h);
    ky4 = RosslerY(*x+h, *y+h*ky3);
    kz4 = RosslerZ(*x+h, *z+h*kz3);
    
    *x = *x+(1.0/6)*h*(kx1+2*kx2+2*kx3+kx4);
    *y = *y+(1.0/6)*h*(ky1+2*ky2+2*ky3+ky4);
    *z = *z+(1.0/6)*h*(kz1+2*kz2+2*kz3+kz4);
    
}

void Runge_Kutta(double *x,double *y,double *z,double h,int temporalDelay)
{
    int i;
    double ky1,kx1,kz1,kx2,ky2,kz2,kx3,ky3,kz3,kx4,ky4,kz4;
    
    kx1 = RosslerX(y[temporalDelay-1], *z);
    ky1 = RosslerY_Delay(*x, y[temporalDelay-1], y[0]);
    kz1 = RosslerZ(*x, *z);
    
    kx2 = RosslerX(y[temporalDelay-1]+0.5*h, *z+0.5*h);
    ky2 = RosslerY_Delay(*x+0.5*h, y[temporalDelay-1]+0.5*h*ky1, y[0]);
    kz2 = RosslerZ(*x+0.5*h, *z+0.5*h*kz1);
    
    kx3 = RosslerX(y[temporalDelay-1]+0.5*h, *z+0.5*h);
    ky3 = RosslerY_Delay(*x+0.5*h, y[temporalDelay-1]+0.5*h*ky2, y[0]);
    kz3 = RosslerZ(*x+0.5*h, *z+0.5*h*kz2);
    
    kx4 = RosslerX(y[temporalDelay-1]+h, *z+h);
    ky4 = RosslerY_Delay(*x+h, y[temporalDelay-1]+h*ky3, y[0]);
    kz4 = RosslerZ(*x+h, *z+h*kz3);
    
    for(i=1;i<temporalDelay;i++){
        y[i-1] = y[i];
    }
    
    *x = *x+(1.0/6)*h*(kx1+2*kx2+2*kx3+kx4);
    y[temporalDelay-1] = y[temporalDelay-1]+(1.0/6)*h*(ky1+2*ky2+2*ky3+ky4);
    *z = *z+(1.0/6)*h*(kz1+2*kz2+2*kz3+kz4);
    
}


double RosslerX(double y,double z)
{
    return -y-z;
}

double RosslerY_Delay(double x,double y,double yDelay)
{
    double F;
    F = k*(yDelay-y);
    if(F>0.2)
        F = 0.2;
    if(F<-0.2)
        F = -0.2;
    
    return x+0.2*y+F;
}

double RosslerY(double x,double y)
{
    return x+0.2*y;
}

double RosslerZ(double x,double z)
{
    return 0.2+z*(x-5.7);
}

double Distance(double x1,double x2,double y1,double y2,double z1,double z2)
{
    return sqrt(pow(y1-y2,2)+pow(x1-x2,2)+pow(z1-z2,2));
}
