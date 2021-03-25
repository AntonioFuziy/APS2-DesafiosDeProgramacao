#include <math.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}


double complex sp[MAX_SIZE/2];
double complex si[MAX_SIZE/2];

double complex tp[MAX_SIZE/2];
double complex ti[MAX_SIZE/2];

int is = 0;
int ip = 0;

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    // for(int i = 0; i < n; i++){
    //     if(i % 2 == 0){
    //         sp[is] = s[i];
    //         is ++;
    //     } else{
    //         si[ip] = s[i];
    //         ip ++;
    //     }
    // }
    for(int k = 0; k < n; k+=2){
        sp[k] = s[k];
    }   

    for(int a = 1; a < n; a+=2){
        si[a] = s[a];
    }
    
    if(n == 1){
        return;
    }
    
    fft(sp, tp, n/2, sign);
    fft(si, ti, n/2, sign);

    for(int j = 0; j < n/2; j++){
        t[j] = tp[j] + ti[j] * cexp(sign * 2 * PI * j * I / n);
        t[j + n/2] = tp[j] - ti[j] * cexp(sign * 2 * PI * j * I / n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}