#include <stdio.h>
#include <math.h>

void c_markov(int *data, double *warm, double *cold, int *size){
    int i=0;
    double warm_count=0;
    double kalt_count=0;    
    double warm_per=0;
    double kalt_per=0;
    for (i=0; i<((*size)-1); ++i)
    {
        if (*(data+i) == 1)
        {
            warm_count = warm_count + 1;
            if (*(data+i) == *(data+i+1))
            {
                warm_per = warm_per+1;
            }
        }
        if (*(data+i) == -1)
        {
            kalt_count = kalt_count + 1;
            if (*(data+i) == *(data+i+1))
            {
                kalt_per = kalt_per+1;
            }
        }
    }

    *warm = warm_per/warm_count;
    *cold = kalt_per/kalt_count;
}

void c_markov_season(int *data, double *warm_w, double *cold_w, double *warm_s, double *cold_s, int *start, int *season_length, int *season_number, int *size){
    int i=*start;
    int j=0;
    int season=1;
    int counter_w=0;
    int counter_s=0;

    double warm_count=0;
    double kalt_count=0;    
    double warm_per=0;
    double kalt_per=0;


    while (i+*(season_length + season)  < *size)
    {   
        if ((*(data+i)!=0) && (*(data+i+*(season_length+season))!=0))
        {
            //printf("%d season %d i %d data\n",season,i,*(data+i) );
            warm_count=0;
            kalt_count=0;    
            warm_per=0;
            kalt_per=0;
            for (j=0; j < (*(season_length+season)); ++j)
            {
                if (*(data+i) == 1)
                {
                    warm_count = warm_count + 1;
                    if (*(data+i) == *(data+i+1))
                    {
                        warm_per = warm_per+1;
                    }
                }
                if (*(data+i) == -1)
                {
                    kalt_count = kalt_count + 1;
                    if (*(data+i) == *(data+i+1))
                    {
                        kalt_per = kalt_per+1;
                    }
                }
                i=i+1;
            }
            if (season==0)
            {
                *(warm_w+counter_w) = warm_per/warm_count;
                *(cold_w+counter_w) = kalt_per/kalt_count;
            }
            if (season==1)
            {
                *(warm_s+counter_s) = warm_per/warm_count;
                *(cold_s+counter_s) = kalt_per/kalt_count;
            }
            
        } else {
            i = i + *(season_length+season);
        }
        if (season==0)
        {
            counter_w=counter_w+1;

        }
        if (season==1)
        {
            counter_s=counter_s+1;
        }

        season=season+1;
        if (season == *season_number){
            season=0;
        }

    }
}

void c_markov_year(int *data, double *warm_y, double *cold_y, int *start, int *interval, int *size){
    int i=*start;
    int j=0;
    int counter=0;

    double warm_count=0;
    double kalt_count=0;    
    double warm_per=0;
    double kalt_per=0;


    while (i+*interval  < *size)
    {   
        if ((*(data+i)!=0) && (*(data+i+*interval)!=0))
        {
            //printf("%d season %d i %d data\n",season,i,*(data+i) );
            warm_count=0;
            kalt_count=0;    
            warm_per=0;
            kalt_per=0;
            for (j=0; j < *interval; ++j)
            {
                if (*(data+i) == 1)
                {
                    warm_count = warm_count + 1;
                    if (*(data+i) == *(data+i+1))
                    {
                        warm_per = warm_per+1;
                    }
                }
                if (*(data+i) == -1)
                {
                    kalt_count = kalt_count + 1;
                    if (*(data+i) == *(data+i+1))
                    {
                        kalt_per = kalt_per+1;
                    }
                }
                i=i+1;
            }
            *(warm_y+counter) = warm_per/warm_count;
            *(cold_y+counter) = kalt_per/kalt_count;

            
        } else {
            i = i + *interval;
        }
        counter=counter+1;
    }
}

void c_markov_run(int *data, double *warm, double *cold, int *size, int *interval){
    int aushol = (*interval-1)/2;
    int i=0;
    int j=0;
    double warm_count=0;
    double kalt_count=0;    
    double warm_per=0;
    double kalt_per=0;
    for (i= aushol; i< *size-aushol; ++i)
    {
        warm_count=0;
        kalt_count=0;    
        warm_per=0;
        kalt_per=0;
        for (j=-aushol; j< (aushol-1); ++j)
        {
            if (*(data+i+j) == 1)
            {
                warm_count = warm_count + 1;
                if (*(data+i+j) == *(data+i+j+1))
                {
                    warm_per = warm_per + 1;
                }
            }
            if (*(data+i+j) == -1)
            {
                kalt_count = kalt_count + 1;
                if (*(data+i+j) == *(data+i+j+1))
                {
                    kalt_per = kalt_per + 1;
                }
            }
        }
        *(warm+i) = warm_per/warm_count;
        *(cold+i) = kalt_per/kalt_count;
    }
}

void c_run_mean2d(double *daten, double *datenex, double *temp_trend, int *size, int *buffer){
    int rows = size[0];
    int cols = size[1];
    int tag = buffer[0];
    int jahr = buffer[1];
    int i = 0;
    int j = 0;

    // daten werden in datenex eingetragen
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
        //printf("%f\n", *(daten+i+j*(rows)));
        *(datenex+i+tag+(j+jahr)*(rows+2*tag)) = *(daten+i+j*(rows));
        }
    }   
    // ende des jahrs wird vor anfang gesetzt
    for (i = 0; i < tag; ++i)
    {
        for (j = 0; j < cols-1; ++j)
        {
        *(datenex+i+(j+1+jahr)*(rows+tag*2)) = *(daten+(i+rows-tag)+j*(rows));
        }
    }
    // anfang vom jahr wird nach ende gesetzt
    for (i = 0; i < tag; ++i)
    {
        for (j = 1; j < cols; ++j)
        {
        *(datenex+(i+rows+tag)+(j-1+jahr)*(rows+2*tag)) = *(daten+(i)+j*(rows));
        }
    }
    // mittelwert wird ausgerechnet
    int ii=0;
    int jj=0;
    double sum=0;
    for (i = tag; i < rows+tag; ++i)
    {
        for (j = jahr; j < cols+jahr; ++j)
        {
            sum=0;
            for (ii = i-tag; ii <= i+tag; ++ii)
            {
                for (jj =j-jahr; jj <= j+jahr; ++jj)
                {   
                    sum=sum+*(datenex+ii+jj*(rows+2*tag));
                }
            }
            *(temp_trend+(i-tag)+(j-jahr)*(rows)) = sum/(float)((tag*2+1)*(jahr*2+1));
            //printf("%f\n", *(trend+(i-tag)+(j-jahr)*(rows)));

        }
    }   

}

void c_mann_kendall(double *daten, double *score, int *start, int *stop){
    int i=0;
    int j=0;
    double S=0;
    int n=*stop-*start;

    double D=n*(n-1)/2;
    for (i = *start; i < (*stop-1); ++i)
    {
        for (j = (i+1); j < (*stop); ++j)
        {
            if ((*(daten+j)- *(daten+i))>0)
            {
                S=S+1;
            }
            if ((*(daten+j)- *(daten+i))<0)
            {
                S=S-1;
            }            
        }
    }
    *score=S/D;  
}