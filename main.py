y_derivada_primeira = 0 # PRÉ INCIALIZAÇÃO DE Y' = W
#EDO PVC MÉTODO DO TIRO
def y_derivada_segunda(x: float, y: float, derivada_primeira): #Y" = C1*(((1+Y')**2)**(1/2)) -> W' = C1*(((1+W)**2)**(1/2))
    c = 0.041**(-1)
    return c*((1+(derivada_primeira)**2)**(1/2))

y0 = 15 # Y(0) = 15
y20 = 10 # Y(20) = 10
#-------------------------------

#EDO PVI RUNGE KUTTA 4ª ORDEM

#CALCULAR DE 0 ATÉ 20 COM PASSO 0.01
x_inicial = 0
x_final = 20
h=0.01 # TAMANHO DO PASSO
#-------------------------------
x_atual = x_inicial # PRÉ INICIALIZAÇÃO DE Xi PARA RUNGE KUTTA
y_atual = y0 # PRÉ INICIALIZAÇÃO DE Yi = Y(0) PARA RUNGE KUTTA
y_prox = 0 # PRÉ INICIALIZAÇÃO DE Yi+1
y_derivada_primeira_atual = 0 # PRÉ INICIALIZAÇÃO DE Wi ONDE W(0) TEM QUE SER UM CHUTE
y_derivada_primeira_prox = 0 #PRÉ INICIALIZAÇÃO DE Wi+1
#-------------------------------

while x_atual < x_final:
    y_K1 = y_derivada_primeira_atual #K1 DE Y
    w_K1 = y_derivada_segunda(0, 0, y_derivada_primeira) #K1 DE W
    y_K2 = y_derivada_primeira_atual + (1/2) * w_K1 # K2 DE Y
    w_K2 = y_derivada_segunda(0 + h / 2, 0 + y_K1 / 2, y_derivada_primeira_atual + w_K1 / 2) # K2 DE W
    y_K3 = y_derivada_primeira_atual + (1/2) * w_K2 # K3 DE Y
    w_K3 = y_derivada_segunda(0 + h / 2, 0 + y_K2 / 2, y_derivada_primeira_atual + w_K2 / 2) # K3 DE W
    y_K4 = y_derivada_primeira_atual + w_K3 # K4 DE Y
    w_K4 = y_derivada_segunda(0 + h, 0 + y_K3, y_derivada_primeira_atual + w_K3) # K4 DE W

    y_prox = y_atual + (1/6) * h * (y_K1 + 2 * (y_K2 + y_K3) + y_K4) # Yi+1 = Yi + 1/6 . h (K1+K4+2(K2+K3))
    y_derivada_primeira_prox = y_derivada_primeira_atual + (1/6) * h * (w_K1 + 2 * (w_K2 + w_K3) + w_K4) # Wi+1 = Wi + 1/6 . h (K1+K4+2(K2+K3))

    x_atual = x_atual + h
    ...