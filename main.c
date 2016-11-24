#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define B        0.1   //[T]
#define d_alpha  90    //[der]
#define d_betta  10    //[deg]
#define d_U      5.0   //[V]
#define enSlit_h 0.001 //[m]
#define enSlit_l 0.04  //[m]
#define exSlit_h 0.001 //[m]
#define exSlit_l 0.04  //[m]

//====    Distributions    ====
float _rand();
float uniform(float, float);
float gauss(float, float);
//====    Additional functions    ====
float min(float*, int);
float max(float*, int);
void  linspace(float, float, float*, int);
void  profile(float*, int, float*, float*, int);
void  writeXY2file(float*, float*, int, int);
//====    MS    ====
void start(float, float*, float*, float*, float*);
void end  (float, float, float, float*);
int  hit  (float*, float);

int main()
{
    srand (time(NULL)); 

    float Radius[] = {0.035, 0.07}; 
    float M[] = {1, 2, 3, 4, 6};    int N_M = sizeof(M)/sizeof(M[0]);
    float Z = 1;
    
    float U[1000], I1[1000], I2[1000]; int N = sizeof(U)/sizeof(U[0]);
    linspace(1, 2700, U, N);

    for (int i = 0; i < N; i++)
    {
        printf("U = %f\n", U[i]);
        int counter1 = 0;
        int counter2 = 0;
        for (int m = 0; m < N_M; m++)
            for (int j = 0; j < 500000; j++)
            {
                float r[3];
                end(U[i], M[m], Z, r);
                if (hit(r, Radius[0]))
                    counter1++;
                if (hit(r, Radius[1]))
                    counter2++;
            }
        I1[i] = counter1;
        I2[i] = counter2;
    }
    
    writeXY2file(U, I1, N, 1);
    writeXY2file(U, I2, N, 2);

    printf("Hello, world!\n");
    return 0;
}

//==================
//====    MS    ====
//==================
void start(float U, float *alpha, float *betta, float *r0, float *Ua)
{
    *alpha = gauss(3.14/2, 3.14*d_alpha/180.0);
    *betta = gauss(0, 3.14*d_betta/180.0);
    r0[0] = enSlit_h*(_rand() - 0.5);
    r0[1] = 0;
    r0[2] = enSlit_l*(_rand() - 0.5);
    *Ua    = gauss(U, d_U);
}

void end(float U, float M, float Z, float *r)
{
    float alpha, betta, r0[3], Ua;
    start(U, &alpha, &betta, r0, &Ua);
    
    float m = 1e-3*M/6e23;
    float q = 1.6e-19*Z;
    float V0 = sqrt(2*q*Ua/m);
    float R = m*V0/(q*B);
    r[0] = r0[0] + 2*R*sin(alpha)*cos(betta);
    r[1] = r0[1];
    r[2] = r0[2] + R*cos(alpha)*(3.14 - 2*betta);
    //printf("r = [%f, %f, %f]\n", r[0], r[1], r[2]);
}

int hit(float *r, float R)
{
    if ((r[0] > (2*R - exSlit_h/2)) && (r[0] < (2*R + exSlit_h/2)))
        if ((r[2] > -exSlit_l/2) && (r[2] < exSlit_l/2))
            return 1;
    return 0;
}

//=============================
//====    Distributions    ====
//=============================
float _rand()
{
    return 0.000001*(0 + rand() % 1000001);
}

float uniform(float start, float stop)
{
    return (start + (stop-start)*_rand());
}

float gauss(float E, float Var)
{
    float x, y, s;
    for(;;)
    {
        x = uniform(-1,1);
        y = uniform(-1,1);
        s = x*x + y*y;
        if (s>0 && s<=1)
            return (E + Var*x*sqrt((-2*log(s))/s));
    }
}

//====================================
//====    Additional functions    ====
//====================================
float min(float *a, int N)
{
    float min = a[0];
    for (int i = 0; i < N; i++)
        if (a[i] < min)
            min = a[i];
    return min;
}

float max(float *a, int N)
{
    float max = a[0];
    for (int i = 0; i < N; i++)
        if (a[i] > max)
            max = a[i];
    return max;
}

void linspace(float _min, float _max, float *s, int Ns)
{
    float step = (_max - _min)/Ns;
    for (int i = 0; i < Ns; i++)
        s[i] = _min + i*step;
}

void profile(float *a, int Na, float *Px, float *Py, int Np)
{
    linspace(min(a, Na), max(a, Na), Px, Np);
    for (int i = 0; i < (Np - 1); i++)
    {
        int counter = 0;
        for (int j = 0; j < Na; j++)
        {
            if ((a[j] >= Px[i]) && (a[j] < Px[i+1]))
                counter ++;
        }
        Py[i] = counter;
    }
    Py[Np-1] = 0;//!!!HARD-CODE!!!
}

void writeXY2file(float *x, float *y, int N, int num)
{
    FILE *fp;
    char name[] = "TestX.txt";
    name[4] = num + '0';
    fp = fopen(name, "w");
    for (int i = 0; i < N; i++)
        fprintf(fp, "%f\t%f\n", x[i], y[i]);
    fclose(fp);
    return;
}
