/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       *
 *                          Ariel Nogueira Kovaljski                          * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define DEBUG 1             /* flag de depuração */

/*====================== Parâmetros a serem ajustados ========================*/


#define nx      200         /* nx: número de células */
#define Delta_t 0.01        /* Δt: passo de tempo (em segundos) */
#define Delta_x 0.2         /* Δx: largura de cada célula (em m) */
#define t_final 15000.0     /* tempo final da simulação (em segundos) */
#define u_bar   0.2e-2      /* ū: velocidade de escoamento (em m/s) */
#define alpha   0.2e-5      /* α: coeficiente de ??? */
#define c_ini   115.5       /* concentração inicial nos volumes da malha */
#define c_inj   7.0e-5      /* concentração de injeção nos vol. da malha */

/*============================================================================*/

void listParameters();
void initializeArray(double arr[], int len, double value);
void calculateQ(double old[], double new[]);
void calculateQDebug(double old[], double new[]);
void printAndSaveResults(double arr[], int len);

int main(void)
{

    double Q_new[nx];   /* array de Q no tempo n+1 */
    double Q_old[nx];   /* array de Q no tempo n */

    if (DEBUG) puts("[INFO] Depuracao ativada\n");

    listParameters();

    /* inicializa os arrays Q para uma concentração c_ini inicial */
    initializeArray(Q_old, nx, c_ini);
    initializeArray(Q_new, nx, c_ini);

    if (DEBUG) calculateQDebug(Q_old, Q_new); else calculateQ(Q_old, Q_new);
    printAndSaveResults(Q_new, nx);

    return 0;
}

void listParameters()
{
    puts("Parametros\n==========");
    puts("Constantes da equacao:");
    printf("Delta_t = %3.2e, Delta_x = %3.2e, u_bar = %3.2e, alpha = %3.2e, "
           "c_inj = %3.2e\n\n", Delta_t, Delta_t, u_bar, alpha, c_inj);
    puts("Constantes da simulacao:");
    printf("nx = %d, t_final = %f, c_ini = %3.2e\n\n",
           nx, t_final, c_ini);
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

        /*************************************************************
         * 
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha 
         */
        new[0] = old[0] - Delta_t/Delta_x 
                 * (
                      u_bar * (2*old[0] - 2*c_inj)
                    - 
                      alpha * (old[1] - 3*old[0] + 2*c_inj) / Delta_x
                   );

        /*************************************************************
         * 
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^     ^     ^     ^     
         *             Para os volumes do centro da malha
         */
        for (x = 1; x < nx - 1; ++x) {
            new[x] = old[x] - Delta_t/Delta_x
                     * (
                          u_bar * (old[x] - old[x-1]) 
                        - 
                          alpha * (old[x-1] - old[x]) / Delta_x
                       );
        }

        /*************************************************************
         * 
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                                             ^
         *             Para o volume da fronteira direita
         *            x possui valor de nx - 1 nesse ponto
         */
        new[x] = old[x] - Delta_t/Delta_x
                 * (
                      u_bar * (old[x] - old[x-1]) 
                    - 
                      alpha * (old[x-1] - old[x]) / Delta_x
                   );

        /* Atualiza array de valores antigos com os novos para a próxima 
           iteração */
        for (x = 0; x < nx; ++x) {
            old[x] = new[x];
        }

        printf("\rCalculando... %3.2f%% concluido (tempo atual: %.2fs)",
               t/t_final * 100, t);
        fflush(stdout);

    } while ( (t += Delta_t) < t_final);

}

/* Função para depuração */
void calculateQDebug(double old[], double new[])
{
    int x, temp;
    double t = 0;

    /* --------------- Depuração --------------- */
    if (DEBUG) {
        puts("Volumes da Malha\n================");
        for (temp = 0; temp < 5; ++temp) {
            /* gera um linspace com 5 elementos */
            printf("%d\t\t", 1 + temp*(nx-1)/(4));
        }
        putchar('\n');
    }
    /* ----------------------------------------- */

    do {
        new[0] = old[0] - Delta_t/Delta_x 
                 * (
                      u_bar * (2*old[0] - 2*c_inj)
                    - 
                      alpha * (old[1] - 3*old[0] + 2*c_inj) / Delta_x
                   );

        /* ------ Depuração ------ */
        if (DEBUG) {
            printf("\r%f\t", new[0]);
        }
        /* ----------------------- */

        /* Para os volumes do centro da malha */
        for (x = 1; x < nx - 1; ++x) {

            new[x] = old[x] - Delta_t/Delta_x
                     * (
                          u_bar * (old[x] - old[x-1]) 
                        - 
                          alpha * (old[x-1] - old[x]) / Delta_x
                       );

            /* -------- Depuração --------- */
            if ( DEBUG && ( x % (nx/4) == 0) ) 
                printf("%f\t", new[x]);
            /* ---------------------------- */
        }

        new[x] = old[x] - Delta_t/Delta_x
                 * (
                      u_bar * (old[x] - old[x-1]) 
                    - 
                      alpha * (old[x-1] - old[x]) / Delta_x
                   );

        /* ---- Depuração ---- */
        if (DEBUG) {
            printf("%f\ttempo atual: %.2fs (%3.2f%% concluido)",
                   new[x], t, t/t_final * 100);
            fflush(stdout);
        }
        /* ------------------- */

        for (x = 0; x < nx; ++x) {
            old[x] = new[x];
        }

    } while ( (t += Delta_t) < t_final);

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
    if ( (results_file = fopen("results.txt", "w")) == NULL) {
        fputs("[ERR] Houve um erro ao escrever o arquivo \"results.txt\"! "
              "Os resultados nao foram salvos.\n", 
              stderr);
        exit(1);
    }

    /* Imprime os resultados no arquivo "results.txt" */
    for (i = 0; i < len; ++i) {
        fprintf(results_file, "%d,%f\n", i + 1, arr[i]);
    }

    fclose(results_file);   /* fecha o arquivo */

    puts("[INFO] Os resultados foram salvos no arquivo \"results.txt\" "
         "no diretorio do programa.");
}