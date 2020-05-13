#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define __SIZEX (1920*10)
const int size_x = __SIZEX;

#define __SIZEY (1920*10)
const int size_y = __SIZEY;

float M[__SIZEX][__SIZEY];

#define __X 0.0
#define __Y 0.0

#define __DX 1.5
#define __DY 1.5
#define __RE_UPPER (__X + __DX)
#define __RE_LOWER (__X - __DX)
#define __IM_UPPER (__Y + __DY)
#define __IM_LOWER (__Y - __DY)

#define __RE_MAP(x,l,u) ((l) + (x)*((u) - (l))/__SIZEX)
#define __IM_MAP(x,l,u) ((l) + (x)*((u) - (l))/__SIZEY)

#define __POWER 2
#define __TH (M_PI * 0.5)
#define __CONST 0.7*(cos(__TH) + I*sin(__TH))
#define __R 1E20
#define __MAX_N 75
const int max_n = __MAX_N;

char* HSLtoRGB(float H, float S, float L, char*r){
    H = H * 360.0;
    float C = (1 - fabs(2*L-1))*S;
    float X = C * (1 - fabs(fmod(H/60.0,2) - 1));
    float m = L - C/2.0;

    float RGB[3] = {0,0,0};
    if (H < 60.0) {
        //RGB = {C, X, 0};
        RGB[0] = C;
        RGB[1] = X;
    }
    else if (H < 120.0){
        RGB[0] = X;
        RGB[1] = C;
    }
    else if (H < 180.0){
        RGB[1] = C;
        RGB[2] = X;
    }
    else if (H < 240.0){
        RGB[1] = X;
        RGB[2] = C;
    }
    else if (H < 300.0){
        RGB[0] = X;
        RGB[2] = C;
    }
    else {
        RGB[0] = C;
        RGB[2] = X;
    }

    for (int i = 0; i < 3; i++)
       r[i] = (int)((RGB[i]+m)*255.0);
    
    return r;
}

typedef struct julia {
    double x,y,dx,dy;
    double complex power;
    double complex c;


    double R;
    int max_n;

    int min,max;

    char* filename;
} julia_t;

int write_image(julia_t* Set){
    Set->min = Set->max_n;
    Set->max = 0;
    clock_t t = clock();
    for (int x = 0; x < size_x; x++)
        for (int y = 0; y < size_y; y++) {
            double complex z = __RE_MAP((double)x,Set->x - Set->dx, Set->x + Set->dx) + __IM_MAP((double)y,Set->y - Set->dy, Set->y + Set->dy)*I;
            M[x][y] = Set->max_n;
            for (int i = 0; i < Set->max_n; i++){
                if (cabs(z) > Set->R) {
                    M[x][y] = i;
                    if (i > Set->max)
                        Set->max = i;
                    if (i < Set->min)
                        Set->min = i;
                    break;
                }
                z = cpow(z,Set->power) + Set->c;
            }
        }
    t = clock() - t;
    printf("%f, Max %d, Min %d\n", (double)t / CLOCKS_PER_SEC, Set->max, Set->min);
    FILE *fp = fopen(Set->filename, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", size_x, size_y);
    for (int x = 0; x < size_x; x++)
        for (int y = 0; y < size_y; y++){
            unsigned char C[3] = {0xff, 0xff, 0xff};
            // Convergent Point
            if (M[x][y] == Set->max_n){
                C[0] = 0;
                C[1] = 0;
                C[2] = 0;
            }
            else {
                double S = fmod(((double)M[x][y] - Set->min)/(Set->max - Set->min), 1.0);
                if (S > 1.0)
                    S = 1.0;
                else if (S < 0)
                    S = 0.0;
                HSLtoRGB(99.0/360.0,1,S,C);
            }
            fwrite(C, 1, 3, fp);
        }
    fclose(fp);
    return 0;
}


int main(){
    julia_t J;
    J.x = 0;
    J.y = 0;
    J.dx = 1.5;
    J.dy = 1.5;
    J.power = 2;
    J.c = 0.5 * cexp(I*M_PI*0.25);
    J.max_n = 1000;
    J.R = 1E300;
    for (int i = 0; i < 1; i++){
        //J.c = 0.5 * cexp(I * (double)(i) / 360.0 * 2.0 * M_PI);
        char filename[50];
        //sprintf(filename, "%03d.ppm", i);
        sprintf(filename, "julia.ppm");
        J.filename = filename;
        write_image(&J);
    }
}



