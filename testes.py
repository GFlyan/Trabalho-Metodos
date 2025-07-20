import numpy as np
import matplotlib.pyplot as plt

# =============================================
# CONSTANTES DO PROBLEMA
# =============================================
C = 0.041  # Constante física do cabo (m^-1)
Y_INICIAL = 15  # Altura em x=0 (m)
Y_FINAL = 10  # Altura em x=20 (m)
X_INICIO = 0  # Ponto inicial (m)
X_FIM = 20  # Ponto final (m)
PASSO = 0.01  # Tamanho do passo para Runge-Kutta
TOLERANCIA = 0.00001  # Tolerância para critério de parada


# =============================================
# MÉTODO DE RUNGE-KUTTA DE 4ª ORDEM
# =============================================
def runge_kutta_4ordem(x, y, w, passo, derivadas):
    """
    Implementa um passo do método de Runge-Kutta de 4ª ordem
    """
    k1_y, k1_w = derivadas(x, y, w)
    k2_y, k2_w = derivadas(x + passo / 2, y + passo * k1_y / 2, w + passo * k1_w / 2)
    k3_y, k3_w = derivadas(x + passo / 2, y + passo * k2_y / 2, w + passo * k2_w / 2)
    k4_y, k4_w = derivadas(x + passo, y + passo * k3_y, w + passo * k3_w)

    y_novo = y + (passo / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y)
    w_novo = w + (passo / 6) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)

    return y_novo, w_novo


# =============================================
# MÉTODO DO TIRO PARA RESOLVER PVC
# =============================================
def metodo_tiro():
    """
    Implementa o método do tiro para resolver o problema de valor de contorno
    """

    def derivadas_cabo(x, y, w):
        """Retorna y' e y'' para o problema do cabo"""
        dydx = w
        dwdx = C * np.sqrt(1 + w ** 2)
        return dydx, dwdx

    # Chutes iniciais para y'(0)
    w_chute1 = 0.0
    w_chute2 = -1.0

    def resolver_pvi(w_chute):
        """Resolve o PVI para um dado chute de y'(0)"""
        x = X_INICIO
        y = Y_INICIAL
        w = w_chute
        pontos = []

        while x <= X_FIM + PASSO:
            pontos.append((x, y))
            y, w = runge_kutta_4ordem(x, y, w, PASSO, derivadas_cabo)
            x += PASSO

        return pontos

    # Primeira tentativa
    pontos1 = resolver_pvi(w_chute1)
    erro1 = pontos1[-1][1] - Y_FINAL

    # Segunda tentativa
    pontos2 = resolver_pvi(w_chute2)
    erro2 = pontos2[-1][1] - Y_FINAL

    # Método da secante para encontrar o y'(0) correto
    while abs(erro2) > TOLERANCIA:
        w_novo = w_chute2 - (erro2 * (w_chute2 - w_chute1)) / (erro2 - erro1)

        w_chute1, w_chute2 = w_chute2, w_novo
        erro1 = erro2

        pontos2 = resolver_pvi(w_chute2)
        erro2 = pontos2[-1][1] - Y_FINAL

    return pontos2, w_chute2


# =============================================
# DIFERENCIAÇÃO NUMÉRICA MANUAL
# =============================================
def diferenciacao_numerica(x_vals, y_vals):
    """
    Calcula derivadas numericamente sem usar bibliotecas avançadas
    """
    n = len(x_vals)
    dy_dx = np.zeros(n)
    d2y_dx2 = np.zeros(n)

    # Primeira derivada (diferenças centradas)
    for i in range(1, n - 1):
        dy_dx[i] = (y_vals[i + 1] - y_vals[i - 1]) / (2 * PASSO)

    # Extremos (diferenças progressivas/regressivas)
    dy_dx[0] = (-3 * y_vals[0] + 4 * y_vals[1] - y_vals[2]) / (2 * PASSO)
    dy_dx[-1] = (3 * y_vals[-1] - 4 * y_vals[-2] + y_vals[-3]) / (2 * PASSO)

    # Segunda derivada (diferenças centradas)
    for i in range(1, n - 1):
        d2y_dx2[i] = (y_vals[i + 1] - 2 * y_vals[i] + y_vals[i - 1]) / (PASSO ** 2)

    # Extremos
    d2y_dx2[0] = (2 * y_vals[0] - 5 * y_vals[1] + 4 * y_vals[2] - y_vals[3]) / (PASSO ** 2)
    d2y_dx2[-1] = (2 * y_vals[-1] - 5 * y_vals[-2] + 4 * y_vals[-3] - y_vals[-4]) / (PASSO ** 2)

    return dy_dx, d2y_dx2


# =============================================
# REGRESSÃO POLINOMIAL MANUAL (4º GRAU)
# =============================================
def regressao_polinomial(x_vals, y_vals, grau=4):
    """
    Implementa regressão polinomial manualmente
    """
    # Monta a matriz de Vandermonde
    A = np.zeros((len(x_vals), grau + 1))
    for i in range(grau + 1):
        A[:, i] = x_vals ** (grau - i)

    # Resolve o sistema por mínimos quadrados
    coeficientes = np.linalg.lstsq(A, y_vals, rcond=None)[0]

    # Função para calcular o polinômio
    def polinomio(x):
        result = 0
        for i in range(grau + 1):
            result += coeficientes[i] * x ** (grau - i)
        return result

    # Calcula derivadas analíticas
    def derivada_polinomio(x):
        result = 0
        for i in range(grau):
            result += (grau - i) * coeficientes[i] * x ** (grau - i - 1)
        return result

    def segunda_derivada_polinomio(x):
        result = 0
        for i in range(grau - 1):
            result += (grau - i) * (grau - i - 1) * coeficientes[i] * x ** (grau - i - 2)
        return result

    return polinomio, derivada_polinomio, segunda_derivada_polinomio


# =============================================
# EXECUÇÃO PRINCIPAL
# =============================================
if __name__ == "__main__":
    # 1. Resolve o problema usando método do tiro
    pontos_solucao, w_inicial = metodo_tiro()
    x_vals = np.array([p[0] for p in pontos_solucao])
    y_vals = np.array([p[1] for p in pontos_solucao])

    print(f"Valor encontrado para y'(0): {w_inicial:.6f}")

    # 2. Plota a solução
    plt.figure(figsize=(12, 6))
    plt.plot(x_vals, y_vals, label='Solução numérica')
    plt.scatter([0, 20], [15, 10], color='red', label='Pontos de fixação')
    plt.title('Forma de um cabo suspenso entre dois pontos')
    plt.xlabel('Distância horizontal (m)')
    plt.ylabel('Altura (m)')
    plt.grid(True)
    plt.legend()

    # 3. Verificação por diferenciação numérica
    dy_dx, d2y_dx2 = diferenciacao_numerica(x_vals, y_vals)
    lado_direito = C * np.sqrt(1 + dy_dx ** 2)
    erro_diff = np.abs(d2y_dx2 - lado_direito)
    print(f"Erro máximo na diferenciação numérica: {np.max(erro_diff):.2e}")

    # 4. Verificação por regressão polinomial
    polinomio, dy_poly, d2y_poly = regressao_polinomial(x_vals, y_vals)

    # Avalia em todos os pontos
    y_poly = np.array([polinomio(x) for x in x_vals])
    dy_poly_vals = np.array([dy_poly(x) for x in x_vals])
    d2y_poly_vals = np.array([d2y_poly(x) for x in x_vals])

    lado_direito_poly = C * np.sqrt(1 + dy_poly_vals ** 2)
    erro_poly = np.abs(d2y_poly_vals - lado_direito_poly)
    print(f"Erro máximo na regressão polinomial: {np.max(erro_poly):.2e}")

    plt.show()