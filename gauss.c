#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

void rowSwap(uint32_t N, uint32_t M, double A[N][M], uint32_t d, uint32_t i) {
    for (int s = d; s < M; s++) {
        double tmp = A[d][s];
        A[d][s] = A[i][s];
        A[i][s] = tmp;
    }
}

void rowSubstract(uint32_t N, uint32_t M, double A[N][M], uint32_t d, uint32_t i) {
    double coefficient = A[i][d]/A[d][d];

    A[i][d] = 0;

    for (int s = d + 1; s < M; s++) {
        A[i][s] -=  coefficient * A[d][s];
    }
}

void triangulate(uint32_t N, uint32_t M, double A[N][M]) {
    uint32_t d = 0;

    while (d < N && d < M) {
        if (A[d][d] == 0) {
            // Hay un cero en la base
            bool swap = false;
            uint32_t i;

            for (i = d + 1; i < N; i++) {
                if (A[i][d] != 0) {
                    // Encontramos una fila más abajo que es distinta de 0
                    swap = true;
                    break;
                }
            }

            if (swap) {
                // Swappear la fila x con la fila i, a partir de la columna y (incluida)
                rowSwap(N, M, A, d, i);
            } else {
                d++;
            }
        } else {
            // Tenemos algo distinto de cero en la base
            for (uint32_t i = d + 1; i < N; i++) {
                rowSubstract(N, M, A, d, i);
            }

            d++;
        }
    }
}

void print(uint32_t N, uint32_t M, double A[N][M]) {
    for (uint32_t i = 0; i < N; i++) {
        for (uint32_t j = 0; j < M; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
}

enum {
    INFINITY,
    NONE,
    SINGLE
} solution;

// Asume que A ya está triangulada
enum solution solveFor(uint32_t N, uint32_t M, double A[N][M], double b[M]) {
    enum solution output = NONE;

    if (N > M) {
        // Single or no solutions
    } else if (N < M) {
        // Infinite or no solutions
    } else {
        // Single or no solution
    }

    return output;
}

int main(int argc, char *argv[]) {
    double A[3][3] = {{5, 6, 2}, {1, 2, 8}, {8, 2, 3}};
    print(3, 3, A);
    triangulate(3, 3, A);
    print(3, 3, A);
}
