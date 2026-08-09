"""
Punto 3 - Modelo de poblacion  dP/dt = P(P-1)(2-P)
Guia de trabajo 1: Campo de pendientes
Juan Jose Laverde Rios

P(t) esta en miles de individuos, t en anios. Se genera el diagrama de
fase (puntos criticos + recta de fase) y las trayectorias P(t) para las
poblaciones iniciales de los literales b), c), d) y e) del enunciado.

Ejecutar con: python punto_3.py
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

P = sp.symbols("P", real=True)
g_expr = P * (P - 1) * (2 - P)
g_np = sp.lambdify(P, g_expr, "numpy")

COLOR = {"estable": "#2e8b57", "inestable": "#c0392b", "semiestable": "#d68910"}


def puntos_criticos_y_clasificacion():
    poly = sp.Poly(sp.expand(g_expr), P)
    raices = sorted(set(round(float(r.evalf()), 6) for r in sp.real_roots(poly)))
    eps = min(np.diff(raices)) / 6 if len(raices) > 1 else 0.2
    resultado = []
    for c in raices:
        izq, der = g_np(c - eps), g_np(c + eps)
        if izq > 0 and der < 0:
            tipo = "estable"
        elif izq < 0 and der > 0:
            tipo = "inestable"
        else:
            tipo = "semiestable"
        resultado.append((c, tipo))
    return resultado


def graficar_fase(criticos, archivo):
    p_min, p_max = -1, 3.5
    ps = np.linspace(p_min, p_max, 1000)
    gs = g_np(ps)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5.5),
                                    gridspec_kw={"width_ratios": [2.3, 1]})

    ax1.plot(ps, gs, color="steelblue", linewidth=1.8)
    ax1.axhline(0, color="black", linewidth=0.8)
    ax1.fill_between(ps, gs, 0, where=(gs > 0), color="#2e8b57", alpha=0.15)
    ax1.fill_between(ps, gs, 0, where=(gs < 0), color="#c0392b", alpha=0.15)
    for c, tipo in criticos:
        ax1.plot(c, 0, "o", color=COLOR[tipo], ms=9, mec="black", mew=0.8, zorder=5)
        ax1.annotate(f"P={c:.0f}\n({tipo})", (c, 0), textcoords="offset points",
                     xytext=(0, -32), ha="center", fontsize=8)
    ax1.set_xlabel("P (miles)")
    ax1.set_ylabel("dP/dt")
    ax1.set_title("dP/dt = P(P-1)(2-P)")
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
        if g_np(medio) > 0:
            ax2.annotate("", xy=(0, medio + largo), xytext=(0, medio - largo),
                         arrowprops=dict(arrowstyle="-|>", color="#2e8b57", lw=2.2))
        else:
            ax2.annotate("", xy=(0, medio - largo), xytext=(0, medio + largo),
                         arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.2))
    for c, tipo in criticos:
        ax2.plot(0, c, "o", color=COLOR[tipo], ms=10, mec="black", mew=0.8, zorder=5)
        ax2.text(0.12, c, f"P={c:.0f}", va="center", fontsize=8)

    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, archivo)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


def graficar_trayectorias(criticos, archivo):
    # literales b) 3000, c) 1500, d) 500, e) 900 individuos -> P en miles
    condiciones = [
        (3.0, "b) P0 = 3000"),
        (1.5, "c) P0 = 1500"),
        (0.5, "d) P0 = 500"),
        (0.9, "e) P0 = 900"),
    ]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for c, tipo in criticos:
        ax.axhline(c, color=COLOR[tipo], linestyle="--", linewidth=1, alpha=0.7)

    t_eval = np.linspace(0, 20, 600)
    for P0, etiqueta in condiciones:
        sol = solve_ivp(lambda t, P: g_np(P), (0, 20), [P0], t_eval=t_eval)
        ax.plot(sol.t, sol.y[0], linewidth=2.2, label=f"{etiqueta}")
        print(f"{etiqueta}: P(20) = {sol.y[0][-1]:.4f} miles")

    ax.set_xlabel("t (anios)")
    ax.set_ylabel("P (miles de individuos)")
    ax.set_title("Trayectorias P(t) para distintas poblaciones iniciales")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, archivo)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


if __name__ == "__main__":
    criticos = puntos_criticos_y_clasificacion()
    print("Puntos criticos:")
    for c, tipo in criticos:
        print(f"  P = {c:.4f}  ->  {tipo}")

    graficar_fase(criticos, "3_fase.png")
    graficar_trayectorias(criticos, "3_trayectorias.png")
