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

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    double complex sp[n/2];
    double complex si[n/2];

    double complex tp[n/2];
    double complex ti[n/2];
    
    int is = 0;
    int ip = 0;

    if(n == 1){
        t[0] = s[0];
        return;
    }    
    
    for(int i = 0; i < n; i++){
        if(i % 2 == 0){
            sp[is] = s[i];
            is ++;
        } else{
            si[ip] = s[i];
            ip ++;
        }
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

// void calculate_matrix_lines(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height, double complex fourier_vector[], int line){
//     for(int i = 0; i < width; i++){
//         matrix[line][i] = fourier_vector[i];
//     }
// }

// void calculate_matrix_columns(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height, double complex fourier_vector[], int column){
//     for(int i = 0; i < height; i++){
//         matrix[i][column] = fourier_vector[i];
//     }
// }

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for(int l = 0; l < height; l++){
        double complex line[width];
        double complex fourier_line[width];

        for(int i = 0; i < width; i++){
            line[i] = matrix[l][i];
        }

        fft_forward(line, fourier_line, width);

        for(int a = 0; a < width; a++){
            matrix[l][a] = fourier_line[a];
        }
    }

    for(int c = 0; c < width; c++){
        double complex column[height];
        double complex fourier_column[height];

        for(int j = 0; j < height; j++){
            column[j] = matrix[j][c];
        }

        fft_forward(column, fourier_column, height);        

        for(int b = 0; b < height; b++){
            matrix[b][c] = fourier_column[b];
        }
    }
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for(int l = 0; l < height; l++){
        double complex line[width];
        double complex fourier_line[width];

        for(int i = 0; i < width; i++){
            line[i] = matrix[l][i];
        }

        fft_inverse(line, fourier_line, width);

        for(int a = 0; a < width; a++){
            matrix[l][a] = fourier_line[a];
        }
    }

    for(int c = 0; c < width; c++){
        double complex column[height];
        double complex fourier_column[height];

        for(int j = 0; j < height; j++){
            column[j] = matrix[j][c];
        }

        fft_inverse(column, fourier_column, height);        

        for(int b = 0; b < height; b++){
            matrix[b][c] = fourier_column[b];
        }
    }
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
