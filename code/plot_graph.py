###############################################################################
#         Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       #
#                           Ariel Nogueira Kovaljski                          #
###############################################################################

import matplotlib.pyplot as plt

x = []
y = []
parameters = {}

# Abre arquivo para leitura
f = open('results.txt', 'r')

for line_number, line in enumerate(f):
    if line_number < 8:
        # Adiciona parâmetros da simulação em um dicionário
        parameters[line.split('=')[0]] = (line.split('=')[1]).split('\n')[0]
    elif line_number == 8:
        # Pula linha separadora
        next(f)
    else:
        # Separa valores na lista `x` e na lista `y`
        x.append( int(line.split(',')[0]) )
        y.append( float( (line.split(',')[1]).split('\n')[0] ) )

print(parameters)

# Configura e exibe o gráfico
plt.grid(True)
plt.suptitle("Quantidade X Célula")
plt.title(rf"$nx = {parameters['nx']}$, " 
          rf"$\Delta t = {parameters['Delta_t']}$, "
          rf"$\Delta_x = {parameters['Delta_x']}$, "
          rf"$t_{{final}} = {parameters['t_final']}$, "
          rf"$\bar{{u}} = {parameters['u_bar']}$, "
          rf"$\alpha = {parameters['alpha']}$, "
          rf"$c_{{ini}} = {parameters['c_ini']}$, "
          rf"$c_{{inj}} = {parameters['c_inj']}$", fontsize=8)
plt.xlabel("célula (i)")
plt.ylabel("Quantidade (Q)")
plt.plot(x,y,'ko-', markerfacecolor='cyan', markeredgecolor='k')
plt.show()