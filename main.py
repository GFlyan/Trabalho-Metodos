#EDO PVI RUNGE KUTTA 4ª ORDEM
derivadaPrimeira = 0
y0 = 15
y2000 = 10 # 2000 PORQUE H = 0.01 = 20/N | N = 20/0.1 = 2000
#-------------------------------

#EDO PVC MÉTODO DO TIRO
C = 0.041
derivadaSegunda = C*(1+(derivadaPrimeira)**2)**(1/2)
#-------------------------------

#CALCULAR DE X0 ATÉ X2000, OU SEJA DE 0 ATÉ 20
x0 = 0
x2000 = 20
h=0.01 # TAMANHO DO PASSO

# RUNGE KUTTA 4ª ORDEM
yatual = y0
xatual = x0
yprox = 0

while xatual < x2000:
    K1 = ... # derivada primeira aplicada
    K2 = ...
    K3 = ...
    K4 = ...

    yprox = yatual + (1/6*K1+2/6*K2+2/6*K3+1/6*K4)*h
#sla