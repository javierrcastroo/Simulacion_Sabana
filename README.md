Aquí tienes el README actualizado completo para copiar y pegar:

---

# Simulación de los Modos Normales de Vibración de una Sábana Elástica mediante Diferencias Finitas

## Tabla de contenidos

1. [Introducción y motivación](#1-introducción-y-motivación)
2. [Descripción física del sistema](#2-descripción-física-del-sistema)
3. [Modelo físico: la ecuación gobernante](#3-modelo-físico-la-ecuación-gobernante)
4. [Condiciones de contorno](#4-condiciones-de-contorno)
5. [Soluciones analíticas: modos normales de vibración](#5-soluciones-analíticas-modos-normales-de-vibración)
6. [Respuesta en frecuencia y resonancia amortiguada](#6-respuesta-en-frecuencia-y-resonancia-amortiguada)
7. [Discretización del dominio espacial](#7-discretización-del-dominio-espacial)
8. [Discretización temporal: esquema en diferencias finitas](#8-discretización-temporal-esquema-en-diferencias-finitas)
9. [Discretización de las condiciones de contorno](#9-discretización-de-las-condiciones-de-contorno)
10. [Condición de estabilidad CFL](#10-condición-de-estabilidad-cfl)
11. [Paralelización con OpenMP](#11-paralelización-con-openmp)
12. [Resumen del algoritmo](#12-resumen-del-algoritmo)
13. [Compilación y uso](#13-compilación-y-uso)

---

## 1. Introducción y motivación

Una sábana extendida y sometida a tensión es un sistema físico continuo cuya dinámica está gobernada por la **ecuación de onda bidimensional**. Cuando una esquina de la sábana es forzada a oscilar armónicamente, la perturbación se propaga por toda la membrana excitando sus modos propios de vibración. Si la frecuencia de forzamiento coincide con alguna frecuencia natural del sistema se produce **resonancia**: la amplitud crece hasta estabilizarse en un valor finito determinado por el amortiguamiento viscoso del aire.

Desde el punto de vista computacional, la integración numérica sobre una malla de $N \times N$ puntos con coste $\mathcal{O}(N^3)$ para un tiempo de simulación fijo constituye un caso de uso ideal para **OpenMP**: el bucle doble de actualización de nodos es perfectamente paralelizable al no haber dependencias de datos entre iteraciones.

---

## 2. Descripción física del sistema

Se considera una sábana delgada de forma cuadrada con lado $L$, dispuesta **horizontalmente** en el plano $XY$. El desplazamiento se produce en la dirección vertical $Z$.

### 2.1 Sistema de referencia

- Dominio espacial: $\Omega = [0, L] \times [0, L]$
- $u(x, y, t)$ [m]: desplazamiento vertical, positivo hacia arriba
- Gravedad despreciada (oscilaciones pequeñas en plano horizontal)

### 2.2 Hipótesis del modelo

- Desplazamientos pequeños: tensión superficial $T$ constante y uniforme durante el movimiento.
- Densidad superficial de masa $\rho$ [kg/m²] uniforme.
- Amortiguamiento viscoso lineal por rozamiento con el aire: fuerza proporcional a $\partial u/\partial t$.

### 2.3 Configuración de los bordes

| Borde | Ecuación | Condición física |
|---|---|---|
| Borde inferior | $y = 0$ | **Fijo** — anclado |
| Borde superior | $y = L$ | **Parcialmente forzado** — esquina $(0, L)$ oscila, resto fijo |
| Borde izquierdo | $x = 0$ | **Libre** — sin fuerza normal externa |
| Borde derecho | $x = L$ | **Libre** — sin fuerza normal externa |

---

## 3. Modelo físico: la ecuación gobernante

### 3.1 Deducción

Sobre un elemento diferencial de sábana $dx \times dy$ actúan dos fuerzas en la dirección $Z$:

**Fuerza elástica de membrana:**

$$dF_\text{elástica} = T\,\nabla^2 u\,dx\,dy$$

**Fuerza de amortiguamiento viscoso (rozamiento con el aire):**

$$dF_\text{aire} = -b\,\dot{u}\,dx\,dy$$

donde $b$ [kg/(m²·s)] es el coeficiente de amortiguamiento viscoso externo.

Aplicando la segunda ley de Newton y dividiendo por $\rho\,dx\,dy$:

$$\frac{\partial^2 u}{\partial t^2} = \frac{T}{\rho}\nabla^2 u - \frac{b}{\rho}\dot{u}$$

### 3.2 Forma canónica

Definiendo $c = \sqrt{T/\rho}$ (velocidad de propagación) y $\beta = b/\rho$ (tasa de amortiguamiento):

$$\boxed{\frac{\partial^2 u}{\partial t^2} = c^2\nabla^2 u - \beta\,\dot{u}}$$

| Parámetro | Definición | Unidades | Significado |
|---|---|---|---|
| $c$ | $\sqrt{T/\rho}$ | m/s | Velocidad de propagación de ondas |
| $\beta$ | $b/\rho$ | 1/s | Tasa de amortiguamiento viscoso externo |

### 3.3 Interpretación física

- $c^2\nabla^2 u$: fuerza restauradora elástica. Gobierna la propagación de ondas.
- $-\beta\,\dot{u}$: amortiguamiento por rozamiento con el aire. Frena el movimiento de todos los puntos uniformemente, acotando la amplitud en resonancia.

---

## 4. Condiciones de contorno

### 4.1 Borde fijo: $y = 0$ (Dirichlet homogéneo)

$$u(x, 0, t) = 0 \quad \forall x \in [0, L],\;\forall t \geq 0$$

### 4.2 Borde superior: $y = L$ (Dirichlet mixto)

El resto del borde está fijo y solo la esquina izquierda $(x=0, y=L)$ oscila armónicamente:

$$u(0, L, t) = A\sin(\omega t), \quad u(x, L, t) = 0 \quad \forall x \in (0, L]$$

Excitar un único punto proyecta energía sobre **todos los modos** simultáneamente, generando patrones de interferencia complejos y físicamente ricos.

### 4.3 Bordes libres: $x = 0$ y $x = L$ (Neumann homogéneo)

$$\frac{\partial u}{\partial x}(0, y, t) = 0, \quad \frac{\partial u}{\partial x}(L, y, t) = 0 \quad \forall y \in [0, L],\;\forall t \geq 0$$

Las ondas se reflejan perfectamente en estos bordes.

### 4.4 Condiciones iniciales

$$u(x, y, 0) = 0, \qquad \frac{\partial u}{\partial t}(x, y, 0) = 0$$

---

## 5. Soluciones analíticas: modos normales de vibración

Para el sistema sin forzamiento ni amortiguamiento, la separación de variables conduce a:

**Formas propias:**

$$\phi_{mn}(x,y) = \cos\!\left(\frac{m\pi x}{L}\right)\sin\!\left(\frac{n\pi y}{L}\right), \quad m = 0,1,2,\ldots;\; n = 1,2,3,\ldots$$

**Frecuencias angulares naturales:**

$$\boxed{\omega_{mn} = \frac{c\pi}{L}\sqrt{m^2 + n^2}}$$

Primeros modos ordenados por frecuencia:

| Modo $(m, n)$ | $\sqrt{m^2+n^2}$ | $f_{mn}/f_{01}$ | Descripción |
|:---:|:---:|:---:|---|
| $(0, 1)$ | $1$ | $1.000$ | Fundamental: uniforme en $x$, una semionda en $y$ |
| $(1, 1)$ | $\sqrt{2}$ | $1.414$ | Una semionda en $x$, una en $y$ |
| $(0, 2)$ | $2$ | $2.000$ | Uniforme en $x$, dos semiondes en $y$ |
| $(1, 2)$ | $\sqrt{5}$ | $2.236$ | Una en $x$, dos en $y$ |
| $(2, 1)$ | $\sqrt{5}$ | $2.236$ | Dos en $x$, una en $y$ |
| $(0, 3)$ | $3$ | $3.000$ | Uniforme en $x$, tres semiondes en $y$ |

Para excitar el modo $(m, n)$ se fija la frecuencia de forzamiento a $\omega = \omega_{mn} = (c\pi/L)\sqrt{m^2+n^2}$.

---

## 6. Respuesta en frecuencia y resonancia amortiguada

### 6.1 Ecuación modal

Proyectando la ecuación gobernante sobre cada modo propio se obtiene un oscilador amortiguado:

$$\ddot{q}_{mn} + 2\zeta_{mn}\omega_{mn}\dot{q}_{mn} + \omega_{mn}^2 q_{mn} = F_{mn}(t)$$

con factor de amortiguamiento modal:

$$\zeta_{mn} = \frac{\beta}{2\omega_{mn}}$$

### 6.2 Amplitud en resonancia

Para forzamiento armónico en resonancia exacta ($\omega = \omega_{mn}$):

$$|q_{mn}|_\text{res} = \frac{F_0/\omega_{mn}^2}{2\zeta_{mn}}$$

La amplitud es **finita** gracias a $\beta > 0$. El sistema evoluciona desde reposo, crece durante el transitorio de duración $\tau = 1/(\zeta_{mn}\omega_{mn})$, y se estabiliza en el valor estacionario. Aumentar $\beta$ (parámetro `B` en el código) amortigua más rápido y reduce la amplitud máxima.

---

## 7. Discretización del dominio espacial

### 7.1 Malla uniforme

$$h = \frac{L}{N-1}, \quad x_i = ih, \quad y_j = jh, \quad i,j = 0,\ldots,N-1$$

Se denota $u_{i,j}^k \approx u(x_i, y_j, k\Delta t)$.

### 7.2 Laplaciano discreto: stencil de 5 puntos

$$\left(\nabla^2 u\right)_{i,j}^k \approx \frac{u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k}{h^2} + \mathcal{O}(h^2)$$

Se define el numerador $\mathcal{L}_{i,j}^k = u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k$.

---

## 8. Discretización temporal: esquema en diferencias finitas

### 8.1 Ensamblaje

Aplicando diferencias centradas de segundo orden para la aceleración y diferencias hacia atrás para la velocidad:

$$\frac{u_{i,j}^{k+1} - 2u_{i,j}^k + u_{i,j}^{k-1}}{\Delta t^2} = \frac{c^2}{h^2}\mathcal{L}_{i,j}^k - \beta\frac{u_{i,j}^k - u_{i,j}^{k-1}}{\Delta t}$$

Despejando $u_{i,j}^{k+1}$:

$$\boxed{u_{i,j}^{k+1} = (2 - \beta\Delta t)\,u_{i,j}^k + (\beta\Delta t - 1)\,u_{i,j}^{k-1} + r^2\,\mathcal{L}_{i,j}^k}$$

donde $r^2 = c^2\Delta t^2/h^2$.

Para $\beta = 0$ se recupera el esquema clásico de la membrana ideal: $u_{i,j}^{k+1} = 2u_{i,j}^k - u_{i,j}^{k-1} + r^2\mathcal{L}_{i,j}^k$.

### 8.2 Propiedades

| Propiedad | Valor |
|---|---|
| Tipo | Explícito (leapfrog con amortiguamiento) |
| Orden temporal | $\mathcal{O}(\Delta t^2)$ elástico, $\mathcal{O}(\Delta t)$ disipativo |
| Orden espacial | $\mathcal{O}(h^2)$ |
| Coste por paso | $\mathcal{O}(N^2)$ |
| Sistema lineal | No |

---

## 9. Discretización de las condiciones de contorno

### 9.1 Borde fijo ($j = 0$)

$$u_{i,0}^k = 0 \quad \forall i,k$$

### 9.2 Borde superior ($j = N-1$)

$$u_{0,N-1}^k = A\sin(\omega\,k\,\Delta t), \quad u_{i,N-1}^k = 0 \quad \forall i > 0$$

### 9.3 Bordes libres: nodos fantasma

La condición $\partial u/\partial x = 0$ se implementa con nodos virtuales fuera del dominio:

$$u_{-1,j}^k = u_{1,j}^k \implies \mathcal{L}_{0,j}^k = 2u_{1,j}^k + u_{0,j+1}^k + u_{0,j-1}^k - 4u_{0,j}^k$$

$$u_{N,j}^k = u_{N-2,j}^k \implies \mathcal{L}_{N-1,j}^k = 2u_{N-2,j}^k + u_{N-1,j+1}^k + u_{N-1,j-1}^k - 4u_{N-1,j}^k$$

---

## 10. Condición de estabilidad CFL

### 10.1 Condición hiperbólica

$$r = \frac{c\,\Delta t}{h} \leq \frac{1}{\sqrt{2}} \approx 0.707$$

### 10.2 Condición de amortiguamiento

$$\beta\,\Delta t < 1$$

### 10.3 Paso temporal de diseño

Se adopta $r = 0.5$ como valor de seguridad:

$$\boxed{\Delta t = \frac{0.5\,h}{c}}$$

El coste total escala como $\mathcal{O}(N^3)$ para tiempo de simulación fijo, lo que justifica la paralelización.

---

## 11. Paralelización con OpenMP

El único bucle que merece paralelización es el bucle doble de actualización de nodos interiores (paso 1), que concentra prácticamente todo el coste computacional.

```c
#pragma omp parallel for schedule(static)
for (int j = 1; j <= N - 2; j++) {
    for (int i = 0; i < N; i++) {
        ...
        u_next[i][j] = (2.0 - beta*dt) * u_curr[i][j]
                     + (beta*dt - 1.0) * u_prev[i][j]
                     + r2 * lap;
    }
}
```

La paralelización es correcta porque `u_next[i][j]` solo lee de `u_curr` y `u_prev` (solo lectura en este paso), sin ninguna dependencia entre iteraciones. No hay carreras de datos.

Los demás bucles (condiciones de contorno, rotación de arrays) son de tamaño $N$ y no se paralelizan: el overhead de lanzar hilos superaría el trabajo.

### Compilación

```bash
# Sin OpenMP
gcc -O2 -o sheet sheet.c -lm

# Con OpenMP
gcc -O2 -fopenmp -o sheet sheet.c -lm

# Controlar número de hilos
OMP_NUM_THREADS=4 ./sheet
```

---

## 12. Resumen del algoritmo

```
PARÁMETROS:
  L, RHO, T_TENS, B, AMP, omega, N, T_SIM

PREPROCESO:
  c    = sqrt(T_TENS / RHO)
  beta = B / RHO
  h    = L / (N-1)
  dt   = 0.5 * h / c
  r2   = (c*dt/h)²
  c1   = 2 - beta*dt
  c2   = beta*dt - 1
  ASSERT beta*dt < 1

INICIALIZACIÓN:
  u_prev = u_curr = 0

BUCLE TEMPORAL k = 1 ... Nstep:

  t = k * dt

  // PASO 1 (paralelo con OpenMP)
  PARA j = 1 hasta N-2:
    PARA i = 0 hasta N-1:
      il = (i==0) ? 1 : i-1
      ir = (i==N-1) ? N-2 : i+1
      lap = u_curr[ir][j] + u_curr[il][j]
          + u_curr[i][j+1] + u_curr[i][j-1]
          - 4*u_curr[i][j]
      u_next[i][j] = c1*u_curr[i][j] + c2*u_prev[i][j] + r2*lap

  // PASO 2: borde fijo
  u_next[i][0] = 0  para todo i

  // PASO 3: esquina forzada
  u_next[0][N-1] = A*sin(omega*t)
  u_next[i][N-1] = 0  para i > 0

  // PASO 4: rotar
  u_prev = u_curr
  u_curr = u_next

FIN BUCLE
```

---

## 13. Compilación y uso

### Requisitos

- GCC con soporte OpenMP
- Python 3 con `numpy`, `matplotlib`, `pillow`

### Pasos

```bash
# 1. Compilar
gcc -O2 -fopenmp -o sheet sheet.c -lm

# 2. Simular (genera frames.bin)
./sheet

# 3. Generar GIF
python make_gif.py
```

### Parámetros ajustables en `sheet.c`

| Parámetro | Efecto |
|---|---|
| `B` | Amortiguamiento. Mayor = más realista, menor amplitud en resonancia |
| `AMP` | Amplitud del forzamiento en la esquina |
| `omega` | Frecuencia de forzamiento. Usar $\omega_{mn} = (c\pi/L)\sqrt{m^2+n^2}$ para resonancia |
| `N` | Resolución de la malla. Mayor = más detalle, más lento |
| `T_SIM` | Tiempo total de simulación |

### Parámetros ajustables en `make_gif.py`

| Parámetro | Efecto |
|---|---|
| `MODE` | `"2d"` (rápido, vista cenital) o `"3d"` (lento, perspectiva) |
| `Z_LIM` | Rango del colormap. Ajustar si la imagen queda saturada |
