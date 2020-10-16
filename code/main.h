/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

/*====================== Parâmetros a serem ajustados ========================*/

#define Lx      100.0       /* Lx: comprimento do domínio (em m) */
#define nx      50          /* nx: número de células */
#define Delta_x (Lx/nx)     /* Δx: largura de cada célula (em m) */
#define u_bar   0.2         /* ū: velocidade de escoamento (em m/s) */
#define alpha   2.0e-4      /* α: coeficiente de difusão */
#define c_ini   1.0         /* concentração inicial nos volumes da malha */
#define c_inj   1.5         /* concentração de injeção nos vol. da malha */
#define t_final 300.0        /* tempo final da simulação (em segundos) */
                            /* Δt: passo de tempo (em segundos) */
#define Delta_t (0.1 * (1/( (2*alpha)/(Delta_x*Delta_x) + u_bar/Delta_x ) ))
/*#define Delta_t 0.1*/

/*============================================================================*/

void listParameters();
void initializeArray(double arr[], int len, double value);
void calculateQ(double old_arr[], double new_arr[]);
void printAndSaveResults(double arr[], int len);