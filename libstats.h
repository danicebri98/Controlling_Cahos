#ifndef _STATS
#define _STATS

#include <math.h>
//FUNCIONES ESTADISTICAS
//MEDIAS
double AverageInt(int *vector,int length)
{
    int i;
    double media;
    media=0;
    for(i=0;i<length;i++)
    {
        media+=vector[i];
    }
    
    return media/length;
}

double AverageDouble(double *vector,int length)
{
    int i;
    double media;
    media=0;
    for(i=0;i<length;i++)
    {
        media+=vector[i];
    }
    media=media/length;
    return media;
}

double AverageFloat(float *vector,int length)
{
    int i;
    float media;
    media=0;
    for(i=0;i<length;i++)
    {
        media+=vector[i];
    }
    media=media/length;
    return media;
}

//VARIANZAS
double VarianceInt(int *vector,int length,double media)
{
    int i;
    double varianza;
    varianza=0;
    for(i=0;i<length;i++)
    {
        varianza+=(vector[i]-media)*(vector[i]-media);
    }
    varianza=varianza/length;
    return varianza;
}

double VarianceDouble(double *vector,int length,double media)
{
    int i;
    double varianza;
    varianza=0;
    for(i=0;i<length;i++)
    {
        varianza+=(vector[i]-media)*(vector[i]-media);
    }
    varianza=varianza/length;
    return varianza;
}

float VarianceFloat(float *vector,int length,double media)
{
    int i;
    float varianza;
    varianza=0;
    for(i=0;i<length;i++)
    {
        varianza+=(vector[i]-media)*(vector[i]-media);
    }
    varianza=varianza/length;
    return varianza;
}

//DESVIACIONES
double DeviationDouble(double varianza,int length)
{
    double desviacion;
    desviacion=sqrt(varianza/length);
    return desviacion;
}

float DeviationFloat(float varianza)
{
    float desviacion;
    desviacion=sqrt(varianza);
    return desviacion;
}

//MOMENTOS N-ÉSIMOS
double SecondMomentFloat(float *vector,int length)
{
    int i;
    float SecMoment;
    
    SecMoment=0;
    
    for(i=0;i<length;i++)
        SecMoment+=vector[i]*vector[i];
    
    return SecMoment/length;
}

//CORRELACIONES
double CorrelationsFloat(float *vector1,float *vector2,int length)
{
    int i;
    float Correlation;
    
    Correlation=0;
    
    for(i=0;i<length;i++)
        Correlation+=vector1[i]*vector2[i];
    
    return Correlation/length;
}

//HISTOGRAMAS
void HistInt(int *vector,int length, int Nbins,int max,int min)
{
    int i,j;
    unsigned int counter;
    int binWidth,Imin,Imax;
    
    binWidth=(1.0*(max-min))/Nbins;
    
    FILE *f;
    f=fopen("Hist09.txt","wt");
    fprintf(f,"Intervalo\tP.Medio\tValor\tError\n");
    
    for(i=0;i<Nbins;i++)
    {
        counter=0;
        Imin=i*binWidth+min;
        Imax=Imin+binWidth;
        for(j=0;j<length;j++)
        {
            if(vector[j]>=Imin && vector[j]<Imax)
            {
                counter++;
            }
        }
        if(counter!=0)
            fprintf(f,"%d\t%d\t%.1f\t%d\t%f\n",Imin,Imax,(Imax+Imin)*0.5,counter,sqrt((1.0*counter/length)*(1-(1.0*counter/length))));
    }
        
    fclose(f);
}

void HistDouble(double *vector,int length, int Nbins,double max,double min)
{
    int i,counter,j;
    double binWidth,Imin,Imax;

    binWidth=(1.0*(max-min))/Nbins;
    
    FILE *f;
    f=fopen("/Users/danixala5/Documents/Master/Métodos de simulación/Assignment10/HistT.txt","wt");
    fprintf(f,"Intervalo\tP.Medio\tValor\tError\n");
    for(i=0;i<Nbins;i++)
    {
        counter=0;
        Imin=i*binWidth+min;
        Imax=Imin+binWidth;
        for(j=0;j<length;j++)
        {
            if(vector[j]>Imin && vector[j]<Imax)
            {
                counter++;
            }
        }
        if(counter!=0)
            fprintf(f,"%.2f\t%.2f\t%.3f\t%f\t%.15f\n",Imin,Imax,(Imax+Imin)*0.5,1.0*counter/(binWidth*length),sqrt(1.0/(binWidth*length)*(1.0*counter/length)*(1-(1.0*counter/length))));
    }
    fclose(f);
}
#endif

