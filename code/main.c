/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "main.h"

int main(void)
{

    double Q_new[nx];   /* array de Q no tempo n+1 */
    double Q_old[nx];   /* array de Q no tempo n */

    puts("\nMNED II - Trabalho 1\n====================");
    puts("por Ariel Nogueira Kovaljski\n");

    listParameters();

    /* inicializa os arrays Q para uma concentração c_ini inicial */
    initializeArray(Q_old, nx, c_ini);
    initializeArray(Q_new, nx, c_ini);

    calculateQ(Q_old, Q_new);
    printAndSaveResults(Q_new, nx);

    return 0;
}

void listParameters()
{
    puts("Parametros\n----------");
    puts("Constantes da equacao:");
    printf("Delta_t = %f, Delta_x = %f, u_bar = %3.2e, alpha = %3.2e, "
           "c_inj = %3.2e\n\n", Delta_t, Delta_t, u_bar, alpha, c_inj);
    puts("Constantes da simulacao:");
    printf("nx = %d, t_final = %f, c_ini = %3.2e\n\n", nx, t_final, c_ini);
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
    int progress = 0, progress_count = 0;
    int progress_incr = (t_final/Delta_t) * 5 / 100;
    double t = 0;

    do {

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha
         */
        new[0] = old[0] - Delta_t/Delta_x * (
                    u_bar * (2*old[0] - 2*c_inj)
                  - alpha * (old[1] - 3*old[0] + 2*c_inj) / Delta_x
                 );

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^     ^     ^     ^
         *             Para os volumes do centro da malha
         */
        for (x = 1; x < nx - 1; ++x) {
            new[x] = old[x] - Delta_t/Delta_x * (
                        u_bar * (old[x] - old[x-1])
                      - alpha * (old[x+1] - 2*old[x] + old[x-1]) / Delta_x
                     );
        }

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                                             ^
         *             Para o volume da fronteira direita
         *            x possui valor de nx - 1 nesse ponto
         */
        new[x] = old[x] - Delta_t/Delta_x * (
                    u_bar * (old[x] - old[x-1])
                  - alpha * (old[x-1] - old[x]) / Delta_x
                 );

        /* incrementa o progresso a cada 5% */
        if (progress_count == progress_incr){
            progress_count = 0;
            ++progress;
            printf("\rCalculando... %d%% concluido", progress * 5);
            fflush(stdout);
        }

        /* Atualiza array de valores antigos com os novos para a próxima
           iteração */
        for (x = 0; x < nx; ++x) {
            old[x] = new[x];
        }

        /* incrementa contador para cada 5% */
        ++progress_count;

    } while ( (t += Delta_t) <= t_final);

}

/* Imprime na tela e salva os resultados num arquivo de saída */
void printAndSaveResults(double arr[], int len)
{
    int i;
    FILE *results_file;     /* Ponteiro para o arquivo de resultados */

    /* Imprime os resultados no console */
    printf("\n\nQ[%d] (tempo final: %.2fs) = [", nx, t_final);
    for (i = 0; i < len - 1; ++i) {
        printf("%f, ", arr[i]);
    }
    printf("%f]\n\n", arr[i]);

    /* Error Handling -- Verifica se é possível criar/escrever o arquivo de
                         resultados */
    if (   (results_file = fopen("./results/results.txt", "w")   ) == NULL
        && (results_file = fopen("./../results/results.txt", "w")) == NULL) {
        fputs("[ERR] Houve um erro ao escrever o arquivo \"results.txt\"! "
              "Os resultados nao foram salvos.\n", stderr);
        exit(1);
    }

    /* Adiciona os resultados no arquivo "results.txt" */
    fprintf(results_file,
            "nx=%d\n"
            "Delta_t=%f\n"
            "Delta_x=%f\n"
            "t_final=%f\n"
            "u_bar=%f\n"
            "alpha=%f\n"
            "c_ini=%f\n"
            "c_inj=%f\n",
            nx, Delta_t, Delta_x, t_final, u_bar, alpha, c_ini, c_inj);
    fputs("********************\n", results_file);

    for (i = 0; i < len; ++i) {
        fprintf(results_file, "%d,%f\n", i + 1, arr[i]);
    }

    fclose(results_file);   /* fecha o arquivo */

    puts("[INFO] Os resultados foram salvos no arquivo \"results.txt\" "
         "no diretorio \"results/\".");
}