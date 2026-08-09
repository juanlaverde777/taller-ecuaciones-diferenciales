"""
Punto 1 - Campo de pendientes
Guia de trabajo 1: Campo de pendientes
Juan Jose Laverde Rios

Para cada literal se dibuja el campo de pendientes de y' = f(x,y), una
familia de curvas solucion (variando la constante de integracion) y la
solucion particular pedida, resaltada.

Ejecutar con: python punto_1.py
Genera las imagenes en la carpeta img/.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

OUT_DIR = os.path.join(os.path.dirname(__file__), "img")
os.makedirs(OUT_DIR, exist_ok=True)


def dibujar_campo(ax, f, xlim, ylim, paso=0.5):
    """Dibuja el campo de direcciones de y' = f(x, y) con flechas de
    longitud uniforme (se normaliza el vector (1, f(x,y)))."""
    xs = np.arange(xlim[0], xlim[1] + paso, paso)
    ys = np.arange(ylim[0], ylim[1] + paso, paso)
    X, Y = np.meshgrid(xs, ys)

    with np.errstate(all="ignore"):
        pendiente = f(X, Y)
    pendiente = np.nan_to_num(pendiente, nan=0.0, posinf=1e6, neginf=-1e6)

    largo = np.hypot(1, pendiente)
    U = 1 / largo
    V = pendiente / largo

    ax.quiver(X, Y, U, V, angles="xy", pivot="mid",
              color="steelblue", alpha=0.55, width=0.0035)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(alpha=0.25)


def guardar(fig, nombre):
    ruta = os.path.join(OUT_DIR, nombre)
    fig.tight_layout()
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    print(f"Guardada {ruta}")


# ---------------------------------------------------------------------
# a) y' = -y - sin(x), y(0) = 1
#    Lineal. Solucion general: y = (cos x - sin x)/2 + C e^(-x)
#    y(0)=1  ->  C = 1/2
# ---------------------------------------------------------------------
def punto_a():
    f = lambda x, y: -y - np.sin(x)
    x = np.linspace(-6, 6, 800)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-6, 6), (-6, 6))

    for C in range(-3, 4):
        y = C * np.exp(-x) + (np.cos(x) - np.sin(x)) / 2
        ax.plot(x, y, color="gray", linewidth=0.9, alpha=0.7)

    y_p = 0.5 * np.exp(-x) + (np.cos(x) - np.sin(x)) / 2
    ax.plot(x, y_p, color="crimson", linewidth=2.4, label="y(0) = 1")
    ax.plot(0, 1, "ko", ms=5)
    ax.set_title("a)  y' = -y - sin(x)")
    ax.legend()
    guardar(fig, "1a.png")


# ---------------------------------------------------------------------
# b) y' = x + y, y(-2) = 2
#    Lineal. Solucion general: y = C e^x - x - 1
#    y(-2)=2  ->  C = e^2
# ---------------------------------------------------------------------
def punto_b():
    f = lambda x, y: x + y
    x = np.linspace(-4, 2.2, 800)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-4, 2.2), (-8, 8))

    for C in [-2, -1, -0.5, 0.5, 1, 2]:
        y = C * np.exp(x) - x - 1
        ax.plot(x, y, color="gray", linewidth=0.9, alpha=0.7)

    y_p = np.exp(2) * np.exp(x) - x - 1
    mascara = np.abs(y_p) < 8.5
    ax.plot(x[mascara], y_p[mascara], color="crimson", linewidth=2.4,
            label="y(-2) = 2")
    ax.plot(-2, 2, "ko", ms=5)
    ax.set_title("b)  y' = x + y")
    ax.legend()
    guardar(fig, "1b.png")


# ---------------------------------------------------------------------
# c) y' = -x^2 + sin(y)
#    No hay condicion inicial en el enunciado -> se elige y(0) = 0.
#    No tiene solucion elemental cerrada, se resuelve numericamente
#    (solve_ivp) hacia adelante y hacia atras en x.
# ---------------------------------------------------------------------
def punto_c():
    f = lambda x, y: -x**2 + np.sin(y)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-3, 3), (-6, 6), paso=0.4)

    for y0 in [-6, -4, -2, 2, 4, 6]:
        sol_f = solve_ivp(f, (0, 3), [y0], dense_output=True, max_step=0.02)
        sol_b = solve_ivp(f, (0, -3), [y0], dense_output=True, max_step=0.02)
        xf = np.linspace(0, sol_f.t[-1], 200)
        xb = np.linspace(sol_b.t[-1], 0, 200)
        ax.plot(xf, sol_f.sol(xf)[0], color="gray", linewidth=0.9, alpha=0.7)
        ax.plot(xb, sol_b.sol(xb)[0], color="gray", linewidth=0.9, alpha=0.7)

    sol_f = solve_ivp(f, (0, 3), [0.0], dense_output=True, max_step=0.02)
    sol_b = solve_ivp(f, (0, -3), [0.0], dense_output=True, max_step=0.02)
    xf = np.linspace(0, sol_f.t[-1], 300)
    xb = np.linspace(sol_b.t[-1], 0, 300)
    ax.plot(xf, sol_f.sol(xf)[0], color="crimson", linewidth=2.4, label="y(0) = 0")
    ax.plot(xb, sol_b.sol(xb)[0], color="crimson", linewidth=2.4)
    ax.plot(0, 0, "ko", ms=5)
    ax.set_title("c)  y' = -x^2 + sin(y)  (sin condicion, se eligio y(0)=0)")
    ax.legend()
    guardar(fig, "1c.png")


# ---------------------------------------------------------------------
# d) (x^2+1) y' + 3xy = 6x   ->  y' = (6x - 3xy)/(x^2+1)
#    Lineal. Factor integrante (x^2+1)^(3/2).
#    Solucion general: y = 2 + C (x^2+1)^(-3/2)
#    Sin condicion inicial en el enunciado -> se elige y(0) = 1  ->  C = -1
# ---------------------------------------------------------------------
def punto_d():
    f = lambda x, y: (6 * x - 3 * x * y) / (x**2 + 1)
    x = np.linspace(-6, 6, 800)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-6, 6), (-4, 8))

    for C in [-6, -3, -1.5, 1.5, 3, 6]:
        y = 2 + C * (x**2 + 1) ** (-1.5)
        ax.plot(x, y, color="gray", linewidth=0.9, alpha=0.7)

    y_p = 2 - (x**2 + 1) ** (-1.5)
    ax.plot(x, y_p, color="crimson", linewidth=2.4, label="y(0) = 1")
    ax.plot(0, 1, "ko", ms=5)
    ax.set_title("d)  (x²+1)y' + 3xy = 6x  (sin condicion, se eligio y(0)=1)")
    ax.legend()
    guardar(fig, "1d.png")


# ---------------------------------------------------------------------
# e) y' = x e^y
#    Separable: e^-y dy = x dx  ->  -e^-y = x^2/2 + C
#    Sin condicion inicial en el enunciado -> se elige y(0) = 0  ->  K=1
#    y = -ln(1 - x^2/2),  valida para |x| < sqrt(2)
# ---------------------------------------------------------------------
def punto_e():
    f = lambda x, y: x * np.exp(y)
    x = np.linspace(-3, 3, 800)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-3, 3), (-4, 4), paso=0.35)

    for K in [0.3, 0.6, 1.5, 2.5, 4]:
        with np.errstate(invalid="ignore"):
            arg = K - x**2 / 2
        valido = arg > 0
        y = -np.log(np.where(valido, arg, np.nan))
        ax.plot(x, y, color="gray", linewidth=0.9, alpha=0.7)

    lim = np.sqrt(2) - 0.01
    xp = np.linspace(-lim, lim, 600)
    y_p = -np.log(1 - xp**2 / 2)
    ax.plot(xp, y_p, color="crimson", linewidth=2.4, label="y(0) = 0")
    ax.plot(0, 0, "ko", ms=5)
    ax.set_title("e)  y' = x e^y  (sin condicion, se eligio y(0)=0)")
    ax.set_ylim(-4, 4)
    ax.legend()
    guardar(fig, "1e.png")


# ---------------------------------------------------------------------
# f) y' = x - y, y(1) = 1
#    Lineal. Solucion general: y = (x - 1) + C e^(-x)
#    y(1)=1  ->  C = e
# ---------------------------------------------------------------------
def punto_f():
    f = lambda x, y: x - y
    x = np.linspace(-4, 6, 800)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    dibujar_campo(ax, f, (-4, 6), (-6, 6))

    for C in [-3, -1.5, -0.5, 0.5, 1.5, 3]:
        y = (x - 1) + C * np.exp(-x)
        ax.plot(x, y, color="gray", linewidth=0.9, alpha=0.7)

    y_p = (x - 1) + np.e * np.exp(-x)
    ax.plot(x, y_p, color="crimson", linewidth=2.4, label="y(1) = 1")
    ax.plot(1, 1, "ko", ms=5)
    ax.set_title("f)  y' = x - y")
    ax.legend()
    guardar(fig, "1f.png")


if __name__ == "__main__":
    punto_a()
    punto_b()
    punto_c()
    punto_d()
    punto_e()
    punto_f()
