"""
Punto 4 - Modelo de poblacion  dP/dt = 3P - 2P^2
Guia de trabajo 1: Campo de pendientes
Juan Jose Laverde Rios

P(t) esta en miles de individuos, t en anios (modelo logistico). Ademas
de las trayectorias de los literales b) y c), se incluye el analisis
simbolico (con sympy) del literal e)/f): nacimientos del 150% por
trimestre (1.5 por individuo por trimestre) y muertes densidad-
dependientes s*P^2 por trimestre, llevado a tasa anual.

Ejecutar con: python punto_4.py
Genera las imagenes en la carpeta img/.
"""

import os

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

OUT_DIR = os.path.join(os.path.dirname(__file__), "img")
os.makedirs(OUT_DIR, exist_ok=True)

P, s = sp.symbols("P s", real=True, positive=False)
g_expr = 3 * P - 2 * P**2
g_np = sp.lambdify(P, g_expr, "numpy")

COLOR = {"estable": "#2e8b57", "inestable": "#c0392b", "semiestable": "#d68910"}


def puntos_criticos_y_clasificacion(expr, var, eps_manual=None):
    poly = sp.Poly(sp.expand(expr), var)
    raices = sorted(set(round(float(r.evalf()), 6) for r in sp.real_roots(poly)))
    f_np = sp.lambdify(var, expr, "numpy")
    eps = eps_manual or (min(np.diff(raices)) / 6 if len(raices) > 1 else 0.2)
    resultado = []
    for c in raices:
        izq, der = f_np(c - eps), f_np(c + eps)
        if izq > 0 and der < 0:
            tipo = "estable"
        elif izq < 0 and der > 0:
            tipo = "inestable"
        else:
            tipo = "semiestable"
        resultado.append((c, tipo))
    return resultado


def graficar_fase(criticos, p_min, p_max, titulo, archivo, expr_np=g_np):
    ps = np.linspace(p_min, p_max, 1000)
    gs = expr_np(ps)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5.5),
                                    gridspec_kw={"width_ratios": [2.3, 1]})

    ax1.plot(ps, gs, color="steelblue", linewidth=1.8)
    ax1.axhline(0, color="black", linewidth=0.8)
    ax1.fill_between(ps, gs, 0, where=(gs > 0), color="#2e8b57", alpha=0.15)
    ax1.fill_between(ps, gs, 0, where=(gs < 0), color="#c0392b", alpha=0.15)
    for c, tipo in criticos:
        ax1.plot(c, 0, "o", color=COLOR[tipo], ms=9, mec="black", mew=0.8, zorder=5)
        ax1.annotate(f"P={c:.3g}\n({tipo})", (c, 0), textcoords="offset points",
                     xytext=(0, -32), ha="center", fontsize=8)
    ax1.set_xlabel("P (miles)")
    ax1.set_ylabel("dP/dt")
    ax1.set_title(titulo)
    ax1.grid(alpha=0.25)

    ax2.set_xlim(-1, 1)
    ax2.set_ylim(p_min, p_max)
    ax2.axvline(0, color="black", linewidth=1.2)
    ax2.set_xticks([])
    ax2.set_title("Recta de fase")

    fronteras = [p_min] + [c for c, _ in criticos] + [p_max]
    for a, b in zip(fronteras[:-1], fronteras[1:]):
        medio = (a + b) / 2
        largo = (b - a) * 0.32
        if expr_np(medio) > 0:
            ax2.annotate("", xy=(0, medio + largo), xytext=(0, medio - largo),
                         arrowprops=dict(arrowstyle="-|>", color="#2e8b57", lw=2.2))
        else:
            ax2.annotate("", xy=(0, medio - largo), xytext=(0, medio + largo),
                         arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.2))
    for c, tipo in criticos:
        ax2.plot(0, c, "o", color=COLOR[tipo], ms=10, mec="black", mew=0.8, zorder=5)
        ax2.text(0.12, c, f"P={c:.3g}", va="center", fontsize=8)

    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, archivo)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


def graficar_trayectorias(criticos, archivo):
    # b) 2000 individuos, c) 100 individuos, d) 1500 individuos (equilibrio)
    condiciones = [
        (2.0, "b) P0 = 2000"),
        (0.1, "c) P0 = 100"),
        (1.5, "d) P0 = 1500 (equilibrio)"),
    ]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for c, tipo in criticos:
        ax.axhline(c, color=COLOR[tipo], linestyle="--", linewidth=1, alpha=0.7)

    t_eval = np.linspace(0, 8, 600)
    for P0, etiqueta in condiciones:
        sol = solve_ivp(lambda t, P: g_np(P), (0, 8), [P0], t_eval=t_eval)
        ax.plot(sol.t, sol.y[0], linewidth=2.2, label=etiqueta)
        print(f"{etiqueta}: P(8) = {sol.y[0][-1]:.4f} miles")

    ax.set_xlabel("t (anios)")
    ax.set_ylabel("P (miles de individuos)")
    ax.set_title("Trayectorias P(t): dP/dt = 3P - 2P^2")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, archivo)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


def literal_e_f():
    """
    Nacimientos: 150% por trimestre -> 1.5*P por trimestre
    Muertes: s*P^2 por trimestre (densidad-dependiente, parametro s)
    Anualizando (4 trimestres/anio): dP/dt = 6P - 4 s P^2
    """
    dPdt_anual = 6 * P - 4 * s * P**2
    print("\nLiteral e): dP/dt =", dPdt_anual, " (anual, con s = tasa de muertes/trimestre)")

    criticos_simbolicos = sp.solve(sp.Eq(dPdt_anual, 0), P)
    print("Puntos criticos simbolicos:", criticos_simbolicos)

    # Grafica ilustrativa con un valor concreto de s (s = 2) para el literal f)
    s_val = 2
    expr_num = dPdt_anual.subs(s, s_val)
    expr_np = sp.lambdify(P, expr_num, "numpy")
    criticos = puntos_criticos_y_clasificacion(expr_num, P)
    graficar_fase(criticos, -0.5, 2, f"e)-f)  dP/dt = 6P - 4sP²  (ilustrado con s = {s_val})",
                  "4ef_fase.png", expr_np=expr_np)

    # Capacidad de carga 3/(2s) en funcion de s
    s_vals = np.linspace(0.2, 5, 200)
    capacidad = 3 / (2 * s_vals)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(s_vals, capacidad, color="steelblue", linewidth=2)
    ax.set_xlabel("s (tasa de muertes por competencia, por trimestre)")
    ax.set_ylabel("Capacidad de carga P* = 3/(2s)  (miles)")
    ax.set_title("Equilibrio estable en funcion de la tasa de muertes s")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, "4ef_capacidad_vs_s.png")
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


if __name__ == "__main__":
    criticos = puntos_criticos_y_clasificacion(g_expr, P)
    print("Puntos criticos:")
    for c, tipo in criticos:
        print(f"  P = {c:.4f}  ->  {tipo}")

    graficar_fase(criticos, -1, 3, "dP/dt = 3P - 2P²", "4_fase.png")
    graficar_trayectorias(criticos, "4_trayectorias.png")
    literal_e_f()
