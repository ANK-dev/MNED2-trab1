/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       *
 *                          Ariel Nogueira Kovaljski                          * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define DEBUG 0                   /* flag de depuração */

/*======================= Parâmetros a serem ajustados =======================*/


#define nx      100         /* nx: número de células */
#define Delta_t 0.01        /* Δt: passo de tempo (em segundos) */
#define Delta_x 0.2         /* Δx: largura de cada célula (em cm) */
#define t_final 6000.0      /* tempo final da simulação (em segundos) */
#define u_bar   0.2         /* ū: velocidade de escoamento (em m/s) */
#define alpha   0.2e-5      /* α: coeficiente de ??? */
#define c_ini   0.5         /* concentração inicial nos volumes da malha */
#define c_inj   0.2         /* concentração de injeção nos vol. da malha */

/*============================================================================*/

/*
 *    +--------   volumes fantasmas   ----------+
 *    |                                         |
 *    v                                         v
 * |--o--|==@==|==@==|==@==|==@==|==@==|==@==|--o--|
 *   i=0   i=1   i=2   i=3   i=4   i=5   i=6   i=7
 *       |<--------------------------------->|
 *                        Lx
 *                volumes principais
 *
 *
 * Equação para a célula da fronteira esquerda no tempo n+1 :
 *                 _                                      _
 *  n+1    n      |  _      n         n      n         n   |
 * Q    = Q  - Δt | 2uc   (t )   - α[Q   - 3Q + 2c   (t )] |
 *  1      1   -- |    inj            2      1    inj      |
 *             Δx |                 ---------------------- |
 *                |_                         Δx           _|
 *
 * Equação para a célula nx no tempo n+1 :
 *                  _                           _
 *  n+1    n       | _  n   n          n     n   |
 * Q    = Q   - Δt | u(Q - Q    ) - α(Q   - Q  ) |
 *  nx     nx   -- |    nx  nx-1       nx-1  nx  |
 *              Δx |                 ----------- |
 *                 |_                     Δx    _|
 *
 * Equação para a célula da fronteira direita no tempo n+1 :
 *                 _                           _
 *  n+1    n      | _  n   n          n     n   |
 * Q    = Q  - Δt | u(Q - Q    ) - α(Q   - Q  ) |
 *  6      6   -- |    6   5          5     6   |
 *             Δx |                 ----------- |
 *                |_                     Δx    _|
 */

void initializeArray(double arr[], int len, double value);
void calculateQ(double old[], double new[]);
void printAndSaveResults(double arr[], int len);

int main(void)
{

    double Q_new[nx];   /* array de Q no tempo n+1 */
    double Q_old[nx];   /* array de Q no tempo n */

    /* inicializa os arrays Q para uma concentração c_ini inicial */
    initializeArray(Q_old, nx, c_ini);
    initializeArray(Q_new, nx, c_ini);

    calculateQ(Q_old, Q_new);
    printAndSaveResults(Q_new, nx);

    return 0;
}

/* Inicializa um array para um valor de entrada */
void initializeArray(double arr[], int len, double value)
{
    int i;
    for (i = 0; i < len; ++i) {
        arr[i] = value;
    }
}

/* Calcula as concentrações na malha ao longo do tempo */
void calculateQ(double old[], double new[])
{
    int x;
    double t = 0;

    do {
        /* --- Depuração --- */
        if (DEBUG) printf("tempo atual: %f\n", t);

        /* Para o volume da fronteira esquerda */
        /* o índice 0 refere-se ao volume nº 1 da malha */
        new[0] = old[0] - Delta_t/Delta_x 
                 * (
                      2*u_bar*c_inj 
                    - 
                      alpha * (old[1] - 3*old[0] + 2*c_inj) / Delta_x
                   );

        /* --- Depuração --- */
        if (DEBUG) printf("%f, ", new[0]);

        for (x = 1; x < nx - 1; ++x) {

            /* Para os volumes do centro da malha */
            new[x] = old[x] - Delta_t/Delta_x
                     * (
                          u_bar * (old[x] - old[x-1]) 
                        - 
                          alpha * (old[x-1] - old[x]) / Delta_x
                       );

            /* --- Depuração --- */
            if (DEBUG) printf("%f, ", new[x]);
        }

        /* Para o volume da fronteira direita */
        /* x possui valor de nx - 1 nesse ponto */
        new[x] = old[x] - Delta_t/Delta_x
                 * (
                      u_bar * (old[x] - old[x-1]) 
                    - 
                      alpha * (old[x-1] - old[x]) / Delta_x
                   );

        if (DEBUG) printf("%f\n\n", new[nx]);

        printf("\rCalculando... %3.2f%% concluido (tempo atual: %.2fs)",
               t/t_final * 100, t);
        fflush(stdout);

    } while ( (t += Delta_t) < t_final);

}

/* Imprime na tela e salva os resultados num arquivo de saída */
void printAndSaveResults(double arr[], int len)
{
    int i;
    FILE *results_file;     /* Ponteiro para o arquivo de resultados */

    /* Imprime os resultados no console */
    printf("\nQ[%d] (tempo final: %.2fs) = [", nx, t_final);

    for (i = 0; i < len - 1; ++i) {
        printf("%f, ", arr[i]);
    }

    printf("%f]\n", arr[i]);

    /* Error Handling -- Verifica se é possível criar/escrever o arquivo de
                         resultados */
    if ( (results_file = fopen("results.txt", "w")) == NULL) {
        fprintf(stderr,
                "[ERR] Houve um erro ao escrever o arquivo \"results.txt\"! "
                "Os resultados nao foram salvos.\n"
               );
        exit(1);
    }

    /* Imprime os resultados no arquivo "results.txt" */
    fprintf(results_file, "\nQ[%d] (tempo final: %.2fs) = [", nx, t_final);

    for (i = 0; i < len - 1; ++i) {
        fprintf(results_file, "%f, ", arr[i]);
    }

    fprintf(results_file, "%f]\n", arr[i]);

    fclose(results_file);   /* fecha o arquivo */

    printf("\n[INFO] Os resultados foram salvos no arquivo \"results.txt\" "
           "no diretório do programa.\n");
}