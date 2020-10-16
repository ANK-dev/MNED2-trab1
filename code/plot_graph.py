###############################################################################
#         Métodos Numéricos para Equações Diferenciais II -- Trabalho 1       #
#                           Ariel Nogueira Kovaljski                          #
###############################################################################

import matplotlib.pyplot as plt

x = []
y = []
parameters = {}

# Abre arquivo para leitura
f = open('./results/results.txt', 'r')

for line_number, line in enumerate(f):
    if line_number < 8:
        # Adiciona parâmetros da simulação em um dicionário
        parameters[line.split('=')[0]] = (line.split('=')[1]).split('\n')[0]
    elif line_number == 8:
        # Pula linha separadora
        pass
    else:
        # Separa valores na lista `x` e na lista `y`
        x.append( int(line.split(',')[0]) )
        y.append( float( (line.split(',')[1]).split('\n')[0] ) )

        # Resolvendo problema da célula
        # x.append( float(line.split(',')[0]) )



# Configura e exibe o gráfico
fig,ax = plt.subplots()
fig.set_size_inches(8, 7)   # Size of the window (1in = 100px)
ax.grid(True)

plt.suptitle("Concentração X Célula")
plt.title(rf"$nx = {parameters['nx']}$, " 
          rf"$\Delta t = {parameters['Delta_t']}$, "
          rf"$\Delta x = {parameters['Delta_x']}$, "
          rf"$t_{{final}} = {parameters['t_final']}$, "
          rf"$\bar{{u}} = {parameters['u_bar']}$, "
          rf"$\alpha = {parameters['alpha']}$, "
          rf"$c_{{ini}} = {parameters['c_ini']}$, "
          rf"$c_{{inj}} = {parameters['c_inj']}$", fontsize=8)
# Rótulo X -- descomentar depois
# plt.xlabel("célula (i)")

# Rótulo X -- resolvendo problema da célula
plt.xlabel("índice célula (i)")
plt.ylabel("concentração (Q)")
plt.plot(x,y,'ko-', markerfacecolor='cyan', markeredgecolor='k')
plt.show()