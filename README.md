# Taller de trabajo 1: Campo de pendientes

**Juan Jose Laverde Rios**
Escuela Colombiana de Ingeniería Julio Garavito
Ecuaciones Diferenciales — 2026

## Cómo lo resolví

Usé GeoGebra Clásico para los campos de direcciones y las curvas solución (comandos `CampoDirecciones` y `ResuelveEDO`), y WolframAlpha (`Roots`, `Reduce`) para hallar y clasificar los puntos críticos en los diagramas de fase, tal como lo indica la guía. Dejo en cada punto los comandos que usé y el análisis correspondiente.

---

## Punto 1 — Campo de pendientes y solución particular

Para cada ecuación dejo el comando del campo de direcciones y el de la curva solución. En los literales que no traían condición inicial elegí un punto arbitrario para poder graficar la familia de soluciones.

### a) y' = -y - sin(x), con y(0) = 1

```
campo_a = CampoDirecciones(-y-sin(x))
ResuelveEDO(-y-sin(x),(0,1))
```

Con el campo se ve que las curvas tienden a "amortiguarse" alrededor de una solución oscilante decreciente en amplitud, coherente con el término `-y` actuando como amortiguamiento sobre el forzamiento `sin(x)`.

### b) y' = x + y, con y(-2) = 2

```
campo_b = CampoDirecciones(x+y)
ResuelveEDO(x+y,(-2,2))
```

Es una lineal de primer orden; la solución que pasa por (-2,2) se aleja rápidamente de la recta y = -x - 1 (que es donde la pendiente se anula), que es la asíntota que se nota en el campo.

### c) y' = -x² + sin(y)

No da condición inicial, así que tomé el punto (0,0).

```
campo_c = CampoDirecciones(-x^2+sin(y))
ResuelveEDO(-x^2+sin(y),(0,0))
```

### d) (x² + 1)y' + 3xy = 6x

Despejando y':

y' = (6x - 3xy) / (x² + 1)

Tomé como condición inicial (0,1).

```
campo_d = CampoDirecciones((6x-3xy)/(x^2+1))
ResuelveEDO((6x-3xy)/(x^2+1),(0,1))
```

### e) y' = x·e^y

Tomé como condición inicial (0,0).

```
campo_e = CampoDirecciones(x*e^y)
ResuelveEDO(x*e^y,(0,0))
```

### f) y' = x - y, con y(1) = 1

```
campo_f = CampoDirecciones(x-y)
ResuelveEDO(x-y,(1,1))
```

La solución que pasa por (1,1) queda prácticamente sobre la recta y = x - 1, que es donde el campo tiene pendiente nula.

---

## Punto 2 — Diagrama de fase y puntos críticos

Para cada ecuación autónoma y' = g(y) hallé las raíces de g con `Roots[...]` y clasifiqué el signo con `Reduce[g(y)>0,y]` / `Reduce[g(y)<0,y]`, igual que en el ejemplo de la guía.

### a) y' = y(3 - y)(y - 2)

```
Roots[c*(3-c)*(c-2)==0,c]
Reduce[c*(3-c)*(c-2)>0,c]
Reduce[c*(3-c)*(c-2)<0,c]
```

Puntos críticos: c = 0, c = 2, c = 3.
g(y) > 0 en (-∞,0) ∪ (2,3); g(y) < 0 en (0,2) ∪ (3,∞).

- y = 0: crece por la izquierda y decrece por la derecha → **asintóticamente estable**.
- y = 2: decrece por la izquierda y crece por la derecha → **inestable**.
- y = 3: crece por la izquierda y decrece por la derecha → **asintóticamente estable**.

### b) y' = y² - y³ = y²(1 - y)

```
Roots[c^2*(1-c)==0,c]
Reduce[c^2*(1-c)>0,c]
Reduce[c^2*(1-c)<0,c]
```

Puntos críticos: c = 0 (raíz doble), c = 1.
g(y) > 0 en (-∞,0) ∪ (0,1); g(y) < 0 en (1,∞).

- y = 0: el signo es positivo a ambos lados (crece por los dos lados) → **semiestable**.
- y = 1: crece por la izquierda, decrece por la derecha → **asintóticamente estable**.

### c) y' = (y + 2)(10 + 3y - y²)

Factorizando 10 + 3y - y² = -(y-5)(y+2), queda g(y) = -(y+2)²(y-5).

```
Roots[(c+2)*(10+3c-c^2)==0,c]
Reduce[(c+2)*(10+3c-c^2)>0,c]
Reduce[(c+2)*(10+3c-c^2)<0,c]
```

Puntos críticos: c = -2 (raíz doble), c = 5.
g(y) > 0 en (-∞,-2) ∪ (-2,5); g(y) < 0 en (5,∞).

- y = -2: crece a ambos lados → **semiestable**.
- y = 5: crece por la izquierda, decrece por la derecha → **asintóticamente estable**.

### d) y' = y⁵ - 4y³ - 5y² = y²(y³ - 4y - 5)

```
Roots[c^5-4c^3-5c^2==0,c]
Reduce[c^5-4c^3-5c^2>0,c]
Reduce[c^5-4c^3-5c^2<0,c]
```

El factor y³ - 4y - 5 solo tiene una raíz real, c ≈ 2.455 (las otras dos son complejas). Los puntos críticos reales son c = 0 (doble) y c ≈ 2.455.
g(y) < 0 en (-∞,0) ∪ (0, 2.455); g(y) > 0 en (2.455,∞).

- y = 0: decrece a ambos lados (a la izquierda se aleja hacia -∞, a la derecha se acerca a 0) → **semiestable**.
- y ≈ 2.455: decrece por la izquierda, crece por la derecha, ambos lados se alejan → **inestable**.

### e) y' = (1 - y)(y - 2)³

```
Roots[(1-c)*(c-2)^3==0,c]
Reduce[(1-c)*(c-2)^3>0,c]
Reduce[(1-c)*(c-2)^3<0,c]
```

Puntos críticos: c = 1, c = 2.
g(y) < 0 en (-∞,1) ∪ (2,∞); g(y) > 0 en (1,2).

- y = 1: decrece por la izquierda, crece por la derecha, ambos lados se alejan de 1 → **inestable**.
- y = 2: crece por la izquierda, decrece por la derecha, ambos lados se acercan a 2 → **asintóticamente estable**.

---

## Punto 3 — Población P(t), dP/dt = P(P-1)(2-P), P en miles

```
Roots[c*(c-1)*(2-c)==0,c]
Reduce[c*(c-1)*(2-c)>0,c]
Reduce[c*(c-1)*(2-c)<0,c]
```

Puntos críticos: P = 0, P = 1, P = 2.
g(P) > 0 en (-∞,0) ∪ (1,2); g(P) < 0 en (0,1) ∪ (2,∞).

### a) Diagrama de fase

- P = 0: **estable** (crece por la izquierda, decrece por la derecha, ambos lados se acercan).
- P = 1: **inestable** (decrece por la izquierda, crece por la derecha, ambos lados se alejan).
- P = 2: **estable** (crece por la izquierda, decrece por la derecha, ambos lados se acercan).

### b) P inicial = 3000 (P₀ = 3)

3 está en la región P > 2, donde g(P) < 0, así que P decrece hacia el equilibrio estable P = 2. A largo plazo la población tiende a **2000 ejemplares**.

### c) P inicial = 1500 (P₀ = 1.5)

1.5 está en (1,2), donde g(P) > 0, así que P crece hacia P = 2. A largo plazo tiende también a **2000 ejemplares**.

### d) P inicial = 500 (P₀ = 0.5)

0.5 está en (0,1), donde g(P) < 0, así que P decrece hacia P = 0. La población **se extingue** con el tiempo.

### e) ¿Puede una población inicial de 900 ejemplares crecer hasta 1100?

P₀ = 0.9 está en (0,1), donde g(P) < 0 (decreciente), y el equilibrio P = 1 es inestable pero actúa como "techo" para esa franja: dentro de (0,1) la población siempre decrece hacia 0. Como 900 < 1000, la trayectoria no puede cruzar el punto crítico inestable P = 1, así que **no**, no puede llegar a 1100; por el contrario tiende a extinguirse (P → 0).

---

## Punto 4 — Población P(t), dP/dt = 3P - 2P², P en miles

g(P) = P(3 - 2P)

```
Roots[3*c-2*c^2==0,c]
Reduce[3*c-2*c^2>0,c]
Reduce[3*c-2*c^2<0,c]
```

Puntos críticos: P = 0, P = 3/2 = 1.5.
g(P) > 0 en (0, 1.5); g(P) < 0 en (-∞,0) ∪ (1.5,∞).

### a) Diagrama de fase

- P = 0: **inestable** (decrece por la izquierda, crece por la derecha, ambos lados se alejan).
- P = 1.5: **asintóticamente estable** (crece por la izquierda, decrece por la derecha, ambos lados se acercan). Es la capacidad de carga del modelo logístico.

### b) P inicial = 2000 (P₀ = 2)

2 está en la región P > 1.5, donde g(P) < 0, así que la población decrece hacia el equilibrio estable. A largo plazo tiende a **1500 ejemplares**.

### c) P inicial = 100 ejemplares (P₀ = 0.1)

0.1 está en (0, 1.5), donde g(P) > 0, así que la población crece hacia el equilibrio estable. A largo plazo tiende también a **1500 ejemplares**.

### d) ¿Qué es correcto afirmar de una población de 1500 ejemplares?

1500 ejemplares corresponde exactamente a P = 1.5, el punto de equilibrio asintóticamente estable del modelo. Es correcto afirmar que esa población **se mantiene constante en el tiempo** (es una solución de equilibrio), y que además es la población hacia la que tienden todas las demás soluciones con P₀ > 0.

### e) Ecuación con nacimientos del 150% por trimestre y muertes s por trimestre

La ecuación original dP/dt = 3P - 2P² está en años, y en ella el término 3P son nacimientos proporcionales a la población y 2P² son muertes por competencia (densidad-dependientes), ambos ya expresados en tasa anual.

Interpreté "tasa de nacimientos del 150% por cada trimestre" como que cada individuo genera en promedio 1.5 individuos nuevos por trimestre, y "mueren s ejemplares en ese mismo periodo" como una muerte densidad-dependiente sP² por trimestre (mismo tipo de término que en el modelo original, pero con tasa s en vez de 2). Como hay 4 trimestres en un año, para llevar la ecuación a tasa anual multiplico ambos términos por 4:

dP/dt = 4(1.5P) - 4(sP²) = 6P - 4sP²

### f) Comportamiento de las soluciones con esa tasa de muertes por trimestre

Con dP/dt = 6P - 4sP² (s > 0), los puntos críticos son P = 0 y P = 6/(4s) = 3/(2s).

```
Roots[6*c-4*s*c^2==0,c]
Reduce[6*c-4*s*c^2>0,c]
```

El comportamiento cualitativo es el mismo del modelo logístico: P = 0 es **inestable** y P = 3/(2s) es **asintóticamente estable**. Es decir, sin importar la población inicial (mientras sea positiva), esta siempre tiende a estabilizarse en 3/(2s) miles de ejemplares. Entre más alta sea la tasa de muertes por competencia s, más baja es la capacidad de carga a la que se estabiliza la población.
