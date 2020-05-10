#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define __SIZEX 3000
const int size_x = __SIZEX;

#define __SIZEY 3000
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

#define __RE_MAP(x) (__RE_LOWER + x*(__RE_UPPER - __RE_LOWER)/__SIZEX)
#define __IM_MAP(x) (__IM_LOWER + x*(__IM_UPPER - __IM_LOWER)/__SIZEY)

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


int main(){
    clock_t t = clock();
    for (int x = 0; x < size_x; x++)
        for (int y = 0; y < size_y; y++) {
            double complex z = __RE_MAP((double)x) + __IM_MAP((double)y)*I;
            M[x][y] = max_n;
            for (int i = 0; i < max_n; i++){
                if (cabs(z) > __R) {
                    M[x][y] = i;
                    //printf("%d\n", i);
                    break;
                }
                z = cpow(z,__POWER) + __CONST;
            }
        }
    t = clock() - t;
    //printf("Computation time: %f\n", (double)t / CLOCKS_PER_SEC);
    FILE *fp = fopen("julia.ppm", "wb");
    fprintf(fp, "P6\n%d %d\n255\n", size_x, size_y);
    for (int x = 0; x < size_x; x++)
        for (int y = 0; y < size_y; y++){
            unsigned char C[3] = {0xff, 0xff, 0xff};
            // Convergent Point
            if (M[x][y] == max_n){
                C[0] = 0;
                C[1] = 0;
                C[2] = 0;
            }
            else {
                double S = (double)M[x][y] / max_n;
                if (S > 1.0)
                    S = 1.0;
                else if (S < 0)
                    S = 0.0;
                HSLtoRGB(S,1,0.5,C);
            }
            fwrite(C, 1, 3, fp);
        }
    fclose(fp);
    return 0;
}






