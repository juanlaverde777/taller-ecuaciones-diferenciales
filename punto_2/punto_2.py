"""
Punto 2 - Diagrama de fase para ecuaciones autonomas
Guia de trabajo 1: Campo de pendientes
Juan Jose Laverde Rios

Para cada ecuacion y' = g(y) se hallan los puntos criticos (raices reales
de g) con sympy, se clasifica su estabilidad evaluando el signo de g a
lado y lado, y se grafica: a la izquierda g(y) contra y (con las regiones
de crecimiento/decrecimiento sombreadas), y a la derecha la recta de fase
con flechas.

Ejecutar con: python punto_2.py
Genera las imagenes en la carpeta img/.
"""

import os

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = os.path.join(os.path.dirname(__file__), "img")
os.makedirs(OUT_DIR, exist_ok=True)

y = sp.symbols("y", real=True)


def puntos_criticos(g_expr):
    poly = sp.Poly(sp.expand(g_expr), y)
    raices = sp.real_roots(poly)
    valores = sorted(set(round(float(r.evalf()), 6) for r in raices))
    return valores


def clasificar(g_np, criticos):
    if len(criticos) > 1:
        gaps = np.diff(criticos)
        eps = max(min(gaps) / 6, 1e-4)
    else:
        eps = 0.2

    resultados = []
    for c in criticos:
        izq = g_np(c - eps)
        der = g_np(c + eps)
        if izq > 0 and der < 0:
            tipo = "estable"
        elif izq < 0 and der > 0:
            tipo = "inestable"
        else:
            tipo = "semiestable"
        resultados.append((c, tipo, izq > 0, der > 0))
    return resultados


COLOR = {"estable": "#2e8b57", "inestable": "#c0392b", "semiestable": "#d68910"}


def graficar(titulo, g_expr, y_min, y_max, archivo):
    g_np = sp.lambdify(y, g_expr, "numpy")
    criticos = puntos_criticos(g_expr)
    clasif = clasificar(g_np, criticos)

    print(f"\n{titulo}")
    for c, tipo, *_ in clasif:
        print(f"  y = {c:.4f}  ->  {tipo}")

    ys = np.linspace(y_min, y_max, 1500)
    gs = g_np(ys)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5.5),
                                    gridspec_kw={"width_ratios": [2.3, 1]})

    ax1.plot(ys, gs, color="steelblue", linewidth=1.8)
    ax1.axhline(0, color="black", linewidth=0.8)
    ax1.fill_between(ys, gs, 0, where=(gs > 0), color="#2e8b57", alpha=0.15)
    ax1.fill_between(ys, gs, 0, where=(gs < 0), color="#c0392b", alpha=0.15)
    for c, tipo, *_ in clasif:
        ax1.plot(c, 0, "o", color=COLOR[tipo], ms=9, mec="black", mew=0.8,
                  zorder=5)
        ax1.annotate(f"y={c:.3g}\n({tipo})", (c, 0), textcoords="offset points",
                     xytext=(0, -32), ha="center", fontsize=8)
    ax1.set_xlabel("y")
    ax1.set_ylabel("g(y)")
    ax1.set_title(titulo)
    ax1.grid(alpha=0.25)

    ax2.set_xlim(-1, 1)
    ax2.set_ylim(y_min, y_max)
    ax2.axvline(0, color="black", linewidth=1.2)
    ax2.set_xticks([])
    ax2.set_title("Recta de fase")

    fronteras = [y_min] + [c for c, *_ in clasif] + [y_max]
    for a, b in zip(fronteras[:-1], fronteras[1:]):
        medio = (a + b) / 2
        signo = g_np(medio)
        largo = (b - a) * 0.32
        if signo > 0:
            ax2.annotate("", xy=(0, medio + largo), xytext=(0, medio - largo),
                         arrowprops=dict(arrowstyle="-|>", color="#2e8b57", lw=2.2))
        else:
            ax2.annotate("", xy=(0, medio - largo), xytext=(0, medio + largo),
                         arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.2))
    for c, tipo, *_ in clasif:
        ax2.plot(0, c, "o", color=COLOR[tipo], ms=10, mec="black", mew=0.8,
                  zorder=5)
        ax2.text(0.12, c, f"y={c:.3g}", va="center", fontsize=8)

    fig.tight_layout()
    ruta = os.path.join(OUT_DIR, archivo)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


if __name__ == "__main__":
    graficar("a)  y' = y(3-y)(y-2)", y * (3 - y) * (y - 2), -2, 5, "2a.png")
    graficar("b)  y' = y^2 - y^3", y**2 - y**3, -2, 3, "2b.png")
    graficar("c)  y' = (y+2)(10+3y-y^2)", (y + 2) * (10 + 3 * y - y**2), -6, 8, "2c.png")
    graficar("d)  y' = y^5 - 4y^3 - 5y^2", y**5 - 4 * y**3 - 5 * y**2, -4, 4, "2d.png")
    graficar("e)  y' = (1-y)(y-2)^3", (1 - y) * (y - 2) ** 3, -1, 4, "2e.png")
