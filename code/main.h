/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

/*====================== Parâmetros a serem ajustados ========================*/

#define Lx      20.0       /* Lx: comprimento do domínio (em m) */
#define nx      50         /* nx: número de células */
#define Delta_x Lx/nx      /* Δx: largura de cada célula (em m) */
#define u_bar   2.0        /* ū: velocidade de escoamento (em m/s) */
#define alpha   2.0e-5     /* α: coeficiente de difusão */
#define c_ini   0.2        /* concentração inicial nos volumes da malha */
#define c_inj   1.2        /* concentração de injeção nos vol. da malha */
#define t_final 10000.0    /* tempo final da simulação (em segundos) */
                           /* Δt: passo de tempo (em segundos) */
#define Delta_t 1/((2*alpha)/(Delta_x*Delta_x) + u_bar/Delta_x)

/*============================================================================*/

void listParameters();
void initializeArray(double arr[], int len, double value);
void calculateQ(double old_arr[], double new_arr[]);
void printAndSaveResults(double arr[], int len);