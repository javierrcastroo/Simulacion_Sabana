# Simulación de los Modos Normales de Vibración de una Sábana Elástica mediante Diferencias Finitas

## Tabla de contenidos

1. [Introducción y motivación](#1-introducción-y-motivación)
2. [Descripción física del sistema](#2-descripción-física-del-sistema)
3. [Ecuación gobernante: la ecuación de onda 2D](#3-ecuación-gobernante-la-ecuación-de-onda-2d)
4. [Condiciones de contorno](#4-condiciones-de-contorno)
5. [Soluciones analíticas: modos normales de vibración](#5-soluciones-analíticas-modos-normales-de-vibración)
6. [Discretización del dominio espacial](#6-discretización-del-dominio-espacial)
7. [Discretización temporal: esquema en diferencias finitas](#7-discretización-temporal-esquema-en-diferencias-finitas)
8. [Discretización de las condiciones de contorno](#8-discretización-de-las-condiciones-de-contorno)
9. [Condición de estabilidad de Courant-Friedrichs-Lewy (CFL)](#9-condición-de-estabilidad-de-courant-friedrichs-lewy-cfl)
10. [Resumen del algoritmo de integración temporal](#10-resumen-del-algoritmo-de-integración-temporal)

---

## 1. Introducción y motivación

Una sábana extendida y sometida a tensión es un sistema físico continuo cuya dinámica está gobernada por una ecuación en derivadas parciales (EDP) hiperbólica: la **ecuación de onda bidimensional**. Cuando uno de sus bordes es forzado a oscilar a una frecuencia determinada, la sábana responde con una combinación de sus modos propios de vibración. Si la frecuencia de forzamiento coincide con alguna de las frecuencias naturales del sistema, se produce **resonancia**, dando lugar a patrones de oscilación estacionarios de gran amplitud conocidos como **modos normales** o **armónicos**.

El estudio de estos modos tiene aplicaciones directas en ingeniería (membranas de altavoces, tímpanos de instrumentos de percusión, velas de embarcaciones), biomecánica (tímpano del oído, tejidos elásticos) y física fundamental (cuerdas de Chladni, figuras de arena sobre placas vibrantes).

Desde el punto de vista computacional, la integración numérica de la ecuación de onda 2D sobre una malla fina constituye un problema de alto coste computacional: en cada paso de tiempo es necesario actualizar el valor de la función en los $N^2$ puntos de la malla, lo que para mallas grandes (e.g. $N = 1000$, con $10^6$ puntos) y miles de pasos temporales supone centenares de millones de operaciones en punto flotante. Este carácter intensivo en cómputo lo convierte en un caso de uso ideal para la **paralelización con OpenMP**.

---

## 2. Descripción física del sistema

Se considera una sábana delgada, flexible e inextensible de forma cuadrada con lado de longitud $L$, sometida a una tensión superficial uniforme $T$ (en N/m). La sábana se modela como una **membrana ideal**, lo que implica las siguientes hipótesis simplificadoras:

- La sábana es perfectamente flexible: no opone resistencia a la curvatura (no hay rigidez a la flexión).
- Los desplazamientos transversales $u$ son **pequeños** en comparación con $L$, de modo que la tensión puede considerarse constante y uniforme en toda la membrana durante el movimiento.
- No existe amortiguamiento (la energía mecánica se conserva).
- La densidad superficial de masa $\rho$ (en kg/m²) es uniforme.

El sistema de referencia se establece de la siguiente manera:

- El plano $xy$ define el plano de reposo de la sábana.
- El dominio espacial es el cuadrado $\Omega = [0, L] \times [0, L]$.
- El desplazamiento transversal fuera del plano en el punto $(x, y)$ y en el instante $t$ se denota por $u(x, y, t)$.

La sábana tiene **cuatro bordes**:

| Borde | Ecuación | Condición física |
|---|---|---|
| Borde inferior | $y = 0$ | **Fijo** (empotrado, nodo) |
| Borde superior | $y = L$ | **Forzado** con movimiento armónico prescrito |
| Borde izquierdo | $x = 0$ | **Libre** (sin restricción transversal) |
| Borde derecho | $x = L$ | **Libre** (sin restricción transversal) |

El borde en $y = 0$ está fijo: no puede desplazarse en ningún instante. El borde en $y = L$, paralelo al anterior, se hace oscilar con una amplitud $A$ y una frecuencia angular $\omega$, actuando como el forzamiento del sistema. Los dos bordes restantes, perpendiculares al eje de excitación, están libres: pueden desplazarse pero no están sometidos a ninguna fuerza normal externa, lo que se traduce en una condición de gradiente nulo en la dirección perpendicular al borde.

---

## 3. Ecuación gobernante: la ecuación de onda 2D

### 3.1 Deducción a partir del principio de Newton

Consideremos un elemento diferencial de sábana de dimensiones $dx \times dy$ y masa $dm = \rho \, dx \, dy$. Las fuerzas transversales ejercidas sobre este elemento por la tensión superficial en las cuatro aristas son:

- En la dirección $x$: la componente transversal neta es $T \, dy \, \left(\frac{\partial^2 u}{\partial x^2}\right) dx$
- En la dirección $y$: la componente transversal neta es $T \, dx \, \left(\frac{\partial^2 u}{\partial y^2}\right) dy$

Aplicando la segunda ley de Newton ($F = ma$) al elemento diferencial:

$$\rho \, dx \, dy \, \frac{\partial^2 u}{\partial t^2} = T \, dy \, \frac{\partial^2 u}{\partial x^2} \, dx + T \, dx \, \frac{\partial^2 u}{\partial y^2} \, dy$$

Dividiendo por $\rho \, dx \, dy$:

$$\frac{\partial^2 u}{\partial t^2} = \frac{T}{\rho} \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)$$

### 3.2 Forma canónica

Definiendo la **velocidad de propagación de ondas** en la membrana:

$$c = \sqrt{\frac{T}{\rho}} \quad \left[\frac{\text{m}}{\text{s}}\right]$$

la ecuación de movimiento adopta su forma canónica:

$$\boxed{\frac{\partial^2 u}{\partial t^2} = c^2 \, \nabla^2 u = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)}$$

Esta es la **ecuación de onda bidimensional** (o ecuación de la membrana), una EDP hiperbólica de segundo orden en el tiempo y en el espacio.

### 3.3 Interpretación física

El operador laplaciano $\nabla^2 u$ mide la curvatura local de la sábana. Cuando la sábana está curvada en un punto, la tensión superficial ejerce una fuerza neta sobre el elemento diferencial en esa zona, que actúa como fuerza restauradora. Esta fuerza es proporcional a $T$ y a la curvatura: a mayor tensión o mayor curvatura, mayor aceleración de retorno. El parámetro $c$ determina con qué rapidez se propagan las perturbaciones a través de la membrana.

---

## 4. Condiciones de contorno

La solución de la ecuación de onda requiere especificar condiciones en todos los bordes del dominio para todo instante $t \geq 0$.

### 4.1 Borde fijo: $y = 0$ (condición de Dirichlet homogénea)

El borde inferior permanece inmóvil en todo momento:

$$u(x, 0, t) = 0 \quad \forall x \in [0, L], \, \forall t \geq 0$$

Esta es una **condición de Dirichlet homogénea**: se prescribe el valor de la función en el contorno. Físicamente corresponde a tener el borde anclado a una estructura rígida.

### 4.2 Borde forzado: $y = L$ (condición de Dirichlet no homogénea)

El borde superior oscila armónicamente con amplitud $A$ y frecuencia angular $\omega$:

$$u(x, L, t) = A \sin(\omega t) \quad \forall x \in [0, L], \, \forall t \geq 0$$

Esta es una **condición de Dirichlet no homogénea**: se prescribe un valor dependiente del tiempo en el contorno. El forzamiento es uniforme a lo largo de todo el borde (no depende de $x$), lo que privilegia la excitación de modos con perfil uniforme en $x$, en particular el modo $m = 0$.

### 4.3 Bordes libres: $x = 0$ y $x = L$ (condición de Neumann homogénea)

Los bordes laterales están libres de fuerza transversal externa. La condición de contorno libre implica que no hay fuerza de tensión en la dirección normal al borde, lo que se traduce en que la derivada normal del desplazamiento es nula:

$$\frac{\partial u}{\partial x}(0, y, t) = 0 \quad \forall y \in [0, L], \, \forall t \geq 0$$

$$\frac{\partial u}{\partial x}(L, y, t) = 0 \quad \forall y \in [0, L], \, \forall t \geq 0$$

Estas son **condiciones de Neumann homogéneas**: se prescribe la derivada normal de la función en el contorno. Físicamente significa que la sábana puede adoptar cualquier altura en los bordes laterales, pero no puede tener pendiente en la dirección $x$ en esos puntos (el borde es siempre localmente horizontal en esa dirección). Matemáticamente, esto equivale a que la sábana se refleja en los bordes $x = 0$ y $x = L$ como si existiesen imágenes simétricas del dominio a ambos lados.

### 4.4 Condiciones iniciales

Para completar el problema de valor inicial y de contorno (PVIC), se necesitan dos condiciones iniciales (dado que la ecuación es de segundo orden en el tiempo):

$$u(x, y, 0) = u_0(x, y) \quad \text{(desplazamiento inicial)}$$

$$\frac{\partial u}{\partial t}(x, y, 0) = v_0(x, y) \quad \text{(velocidad inicial)}$$

En nuestra simulación tomamos $u_0 = 0$ y $v_0 = 0$: la sábana parte del reposo en su posición de equilibrio.

---

## 5. Soluciones analíticas: modos normales de vibración

### 5.1 Separación de variables

Para encontrar las soluciones propias del sistema (sin forzamiento, solo con las condiciones de contorno homogéneas en $y=0$ y $y=L$, y las condiciones de Neumann en $x$), se buscan soluciones de la forma:

$$u(x, y, t) = X(x) \cdot Y(y) \cdot T(t)$$

Sustituyendo en la ecuación de onda y separando variables:

$$\frac{T''}{c^2 T} = \frac{X''}{X} + \frac{Y''}{Y} = -\lambda^2$$

Esto conduce a tres ecuaciones diferenciales ordinarias (EDOs) desacopladas:

$$T'' + c^2 \lambda^2 T = 0$$
$$X'' + k_x^2 X = 0, \quad Y'' + k_y^2 Y = 0$$

con la relación de dispersión: $\lambda^2 = k_x^2 + k_y^2$.

### 5.2 Solución para $X(x)$: condiciones de Neumann

Con las condiciones $X'(0) = 0$ y $X'(L) = 0$, la solución general $X(x) = A\cos(k_x x) + B\sin(k_x x)$ impone $B = 0$ (de $X'(0) = 0$) y $k_x = m\pi/L$ con $m = 0, 1, 2, \ldots$ (de $X'(L) = 0$). Por tanto:

$$X_m(x) = \cos\left(\frac{m\pi x}{L}\right), \quad m = 0, 1, 2, \ldots$$

### 5.3 Solución para $Y(y)$: condiciones de Dirichlet

Con las condiciones $Y(0) = 0$ e $Y(L) = 0$, la solución impone $k_y = n\pi/L$ con $n = 1, 2, 3, \ldots$ y:

$$Y_n(y) = \sin\left(\frac{n\pi y}{L}\right), \quad n = 1, 2, 3, \ldots$$

### 5.4 Modos normales y frecuencias propias

Las soluciones propias del sistema son:

$$u_{mn}(x, y, t) = \cos\left(\frac{m\pi x}{L}\right) \sin\left(\frac{n\pi y}{L}\right) \cos(\omega_{mn} t + \phi_{mn})$$

con **frecuencias angulares propias**:

$$\boxed{\omega_{mn} = c\pi \sqrt{\left(\frac{m}{L}\right)^2 + \left(\frac{n}{L}\right)^2} = \frac{c\pi}{L}\sqrt{m^2 + n^2}}$$

y **frecuencias propias** (en Hz):

$$f_{mn} = \frac{\omega_{mn}}{2\pi} = \frac{c}{2L}\sqrt{m^2 + n^2}$$

### 5.5 Significado de los índices $(m, n)$

- El índice $m$ cuenta el número de **semiondes en la dirección $x$** (nodos verticales internos).
- El índice $n$ cuenta el número de **semiondes en la dirección $y$** (nodos horizontales internos, sin contar los bordes fijo y libre).

Los primeros modos ordenados por frecuencia creciente son:

| Modo $(m, n)$ | $\sqrt{m^2+n^2}$ | Frecuencia relativa $f_{mn}/f_{01}$ | Descripción |
|:---:|:---:|:---:|---|
| $(0, 1)$ | $1$ | $1.000$ | Modo fundamental: oscilación uniforme en $x$, una semionda en $y$ |
| $(1, 1)$ | $\sqrt{2}$ | $1.414$ | Una semionda en $x$, una en $y$ |
| $(0, 2)$ | $2$ | $2.000$ | Dos semiondes en $y$, uniforme en $x$ |
| $(1, 2)$ | $\sqrt{5}$ | $2.236$ | Una en $x$, dos en $y$ |
| $(2, 1)$ | $\sqrt{5}$ | $2.236$ | Dos en $x$, una en $y$ |
| $(2, 2)$ | $2\sqrt{2}$ | $2.828$ | Dos en $x$, dos en $y$ |
| $(0, 3)$ | $3$ | $3.000$ | Tres semiondes en $y$, uniforme en $x$ |

Cuando la frecuencia de forzamiento $\omega$ coincide con algún $\omega_{mn}$, el sistema entra en **resonancia** y la amplitud de oscilación crece de forma sostenida (indefinidamente en ausencia de amortiguamiento).

---

## 6. Discretización del dominio espacial

### 6.1 Malla uniforme

El dominio continuo $\Omega = [0, L] \times [0, L]$ se discretiza mediante una **malla uniforme** de $N \times N$ puntos:

$$x_i = i \cdot h, \quad i = 0, 1, \ldots, N-1$$
$$y_j = j \cdot h, \quad j = 0, 1, \ldots, N-1$$

donde el **paso espacial** es:

$$h = \frac{L}{N - 1}$$

El valor de la solución en el nodo $(i, j)$ y en el instante $k$ se denota:

$$u_{i,j}^k \approx u(x_i, y_j, t_k)$$

### 6.2 Aproximación del laplaciano: diferencias centradas de segundo orden

Se aproxima cada derivada segunda mediante la fórmula de **diferencias finitas centradas**:

$$\frac{\partial^2 u}{\partial x^2}\bigg|_{i,j}^k \approx \frac{u_{i+1,j}^k - 2u_{i,j}^k + u_{i-1,j}^k}{h^2} + \mathcal{O}(h^2)$$

$$\frac{\partial^2 u}{\partial y^2}\bigg|_{i,j}^k \approx \frac{u_{i,j+1}^k - 2u_{i,j}^k + u_{i,j-1}^k}{h^2} + \mathcal{O}(h^2)$$

El error de truncación es de orden $\mathcal{O}(h^2)$: al reducir el paso espacial a la mitad, el error espacial se reduce en un factor de 4.

El laplaciano discreto completo es:

$$\nabla^2 u_{i,j}^k \approx \frac{u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k}{h^2}$$

Este operador es el conocido como **stencil de 5 puntos** (o estrella de 5 puntos), que utiliza el nodo central y sus cuatro vecinos inmediatos (norte, sur, este, oeste).

### 6.3 Validez de la aproximación

La fórmula de diferencias centradas para la segunda derivada se obtiene combinando los desarrollos de Taylor de $u(x+h)$ y $u(x-h)$:

$$u(x+h) = u(x) + h u'(x) + \frac{h^2}{2} u''(x) + \frac{h^3}{6} u'''(x) + \frac{h^4}{24} u^{(4)}(x) + \ldots$$

$$u(x-h) = u(x) - h u'(x) + \frac{h^2}{2} u''(x) - \frac{h^3}{6} u'''(x) + \frac{h^4}{24} u^{(4)}(x) + \ldots$$

Sumando ambas expresiones:

$$u(x+h) + u(x-h) = 2u(x) + h^2 u''(x) + \frac{h^4}{12} u^{(4)}(x) + \ldots$$

Despejando $u''(x)$:

$$u''(x) = \frac{u(x+h) - 2u(x) + u(x-h)}{h^2} - \frac{h^2}{12} u^{(4)}(x) + \ldots$$

El error dominante es $\mathcal{O}(h^2)$, lo que confirma la consistencia de segundo orden del esquema.

---

## 7. Discretización temporal: esquema en diferencias finitas

### 7.1 Aproximación de la derivada temporal de segundo orden

De forma análoga al caso espacial, se aproxima la derivada segunda en el tiempo mediante diferencias centradas:

$$\frac{\partial^2 u}{\partial t^2}\bigg|_{i,j}^k \approx \frac{u_{i,j}^{k+1} - 2u_{i,j}^k + u_{i,j}^{k-1}}{\Delta t^2} + \mathcal{O}(\Delta t^2)$$

donde $\Delta t$ es el **paso temporal** y $u_{i,j}^{k-1}$, $u_{i,j}^k$, $u_{i,j}^{k+1}$ representan los valores en los instantes anteriores, actual y siguiente, respectivamente.

### 7.2 Esquema explícito de integración temporal (FTCS)

Sustituyendo las aproximaciones de la derivada temporal y del laplaciano en la ecuación de onda:

$$\frac{u_{i,j}^{k+1} - 2u_{i,j}^k + u_{i,j}^{k-1}}{\Delta t^2} = c^2 \cdot \frac{u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k}{h^2}$$

Despejando $u_{i,j}^{k+1}$, que es el valor desconocido en el instante siguiente:

$$\boxed{u_{i,j}^{k+1} = 2u_{i,j}^k - u_{i,j}^{k-1} + r^2 \left( u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k \right)}$$

donde $r$ es el **número de Courant** (o número CFL):

$$r = \frac{c \, \Delta t}{h}$$

### 7.3 Propiedades del esquema

Este esquema se denomina **explícito** porque el valor en el instante $k+1$ se calcula directamente a partir de los valores conocidos en los instantes $k$ y $k-1$, sin necesidad de resolver ningún sistema de ecuaciones. Las ventajas e inconvenientes son:

| Característica | Valor |
|---|---|
| Tipo | Explícito (leapfrog) |
| Orden en el tiempo | $\mathcal{O}(\Delta t^2)$ |
| Orden en el espacio | $\mathcal{O}(h^2)$ |
| Coste por paso | $\mathcal{O}(N^2)$ operaciones |
| Necesidad de resolver sistema lineal | No |
| Condición de estabilidad | Sí (véase §9) |

El coste computacional por paso de tiempo es $\mathcal{O}(N^2)$ operaciones, ya que para cada uno de los $N^2$ nodos de la malla se realizan un número fijo de operaciones aritméticas (sumas y multiplicaciones). Este bucle doble sobre los nodos interiores es el núcleo computacional del programa y el objetivo principal de la paralelización.

### 7.4 Inicialización: el primer paso temporal

El esquema necesita los valores en dos instantes anteriores para avanzar: $u^{k-1}$ y $u^k$. Para el primer paso ($k=0 \to k=1$) no se dispone de $u^{-1}$. Se resuelve utilizando la condición inicial de velocidad $v_0 = \partial u / \partial t |_{t=0} = 0$:

$$\frac{u_{i,j}^1 - u_{i,j}^{-1}}{2\Delta t} = 0 \implies u_{i,j}^{-1} = u_{i,j}^1$$

Sustituyendo en la fórmula de actualización para $k=0$:

$$u_{i,j}^1 = u_{i,j}^0 + \frac{r^2}{2} \left( u_{i+1,j}^0 + u_{i-1,j}^0 + u_{i,j+1}^0 + u_{i,j-1}^0 - 4u_{i,j}^0 \right)$$

Con $u_{i,j}^0 = 0$, resulta simplemente $u_{i,j}^1 = 0$ para los nodos interiores (excepto los bordes con condiciones de Dirichlet no homogéneas).

---

## 8. Discretización de las condiciones de contorno

### 8.1 Borde fijo ($j = 0$): Dirichlet homogéneo

La implementación es directa: se impone explícitamente:

$$u_{i,0}^k = 0 \quad \forall i, k$$

Estos nodos nunca se actualizan con la fórmula general; su valor siempre es cero.

### 8.2 Borde forzado ($j = N-1$): Dirichlet no homogéneo

De forma análoga, en cada paso de tiempo $k$ se impone:

$$u_{i,N-1}^k = A \sin(\omega \, t_k) = A \sin(\omega \, k \, \Delta t) \quad \forall i$$

Esta condición actúa como fuente de energía del sistema.

### 8.3 Bordes libres ($i = 0$ e $i = N-1$): Neumann homogéneo

La condición $\partial u / \partial x = 0$ en los bordes laterales no puede aplicarse directamente a los nodos del borde con el stencil de 5 puntos, ya que el nodo $u_{-1,j}^k$ (o $u_{N,j}^k$) no existe en la malla.

Se recurre al método de los **nodos fantasma** (*ghost nodes*): se introducen formalmente nodos virtuales fuera del dominio, $u_{-1,j}^k$ y $u_{N,j}^k$, y se determina su valor a partir de la condición de Neumann mediante diferencias centradas:

$$\frac{u_{1,j}^k - u_{-1,j}^k}{2h} = 0 \implies u_{-1,j}^k = u_{1,j}^k$$

$$\frac{u_{N,j}^k - u_{N-2,j}^k}{2h} = 0 \implies u_{N,j}^k = u_{N-2,j}^k$$

Esto equivale a **reflejar simétricamente** la malla en los bordes laterales. Sustituyendo en la fórmula de actualización para los nodos del borde $i=0$:

$$u_{0,j}^{k+1} = 2u_{0,j}^k - u_{0,j}^{k-1} + r^2 \left( 2u_{1,j}^k + u_{0,j+1}^k + u_{0,j-1}^k - 4u_{0,j}^k \right)$$

Y de forma análoga para $i = N-1$:

$$u_{N-1,j}^{k+1} = 2u_{N-1,j}^k - u_{N-1,j}^{k-1} + r^2 \left( 2u_{N-2,j}^k + u_{N-1,j+1}^k + u_{N-1,j-1}^k - 4u_{N-1,j}^k \right)$$

---

## 9. Condición de estabilidad de Courant-Friedrichs-Lewy (CFL)

### 9.1 Análisis de estabilidad de Von Neumann

El análisis de estabilidad de Von Neumann consiste en estudiar la evolución de un modo de Fourier arbitrario $u_{i,j}^k = \xi^k e^{i(k_x x_i + k_y y_j)}$ bajo la acción del esquema de diferencias finitas. Sustituyendo en la ecuación discretizada y simplificando, se obtiene una ecuación para el **factor de amplificación** $\xi$.

Para que el esquema sea estable (las perturbaciones no crezcan sin límite), es necesario que $|\xi| \leq 1$. El análisis detallado conduce a la condición:

$$r^2 \left(4\sin^2\frac{k_x h}{2} + 4\sin^2\frac{k_y h}{2}\right) \leq 4$$

El caso más restrictivo se da cuando ambas frecuencias espaciales son máximas ($k_x h = k_y h = \pi$), lo que da:

$$r^2 \cdot 4 \leq 4 \implies r \leq 1$$

### 9.2 Condición CFL en 2D

Sin embargo, en 2D la condición más restrictiva proviene del caso en que la perturbación se propaga en dirección diagonal ($k_x = k_y$), lo que lleva a:

$$\boxed{r = \frac{c \, \Delta t}{h} \leq \frac{1}{\sqrt{2}} \approx 0.707}$$

Esta es la **condición de Courant-Friedrichs-Lewy (CFL)** para la ecuación de onda 2D con el esquema explícito de diferencias centradas. Si no se cumple, el error numérico crece exponencialmente con el tiempo y la simulación diverge.

### 9.3 Implicaciones prácticas

La condición CFL establece una relación entre el paso espacial $h$ y el paso temporal máximo permitido:

$$\Delta t \leq \frac{h}{c\sqrt{2}}$$

En la práctica se usa un valor de seguridad $r = 0.5$, lo que da:

$$\Delta t = \frac{0.5 \, h}{c}$$

Esto implica que al refinar la malla espacial por un factor de 2 (doblar $N$, reducir $h$ a la mitad), también hay que reducir $\Delta t$ a la mitad, lo que **cuadruplica el número total de operaciones**: la malla pasa de $N^2$ a $4N^2$ puntos, y el número de pasos temporales se dobla. En total, el coste computacional escala como $\mathcal{O}(N^3)$ para un tiempo de simulación fijo.

Este escalado cúbico es la principal justificación del uso de **paralelización**: para mallas de $N = 2000$ y tiempos de simulación largos, el tiempo de cómputo en serie puede ser de horas, mientras que con paralelización OpenMP se puede reducir a minutos.

---

## 10. Resumen del algoritmo de integración temporal

A continuación se presenta el pseudocódigo completo del algoritmo de simulación, que sintetiza todos los elementos anteriores:

```
PARÁMETROS:
  L     = longitud del lado de la sábana [m]
  N     = número de nodos por lado
  c     = velocidad de onda [m/s]
  A     = amplitud del forzamiento [m]
  omega = frecuencia angular del forzamiento [rad/s]
  T_sim = tiempo total de simulación [s]

PREPROCESO:
  h     = L / (N - 1)             // Paso espacial
  dt    = 0.5 * h / c             // Paso temporal (CFL con r = 0.5)
  r     = c * dt / h              // Número de Courant (= 0.5)
  r2    = r * r                   // r² (precalculado)
  Nstep = ceil(T_sim / dt)        // Número total de pasos

INICIALIZACIÓN:
  u_prev[i][j] = 0  para todo (i, j)   // u en t = -dt (simetría)
  u_curr[i][j] = 0  para todo (i, j)   // u en t = 0

BUCLE TEMPORAL (k = 1, 2, ..., Nstep):

  t_k = k * dt

  // PASO 1: Calcular u_next en nodos interiores
  PARA j = 1 hasta N-2:
    PARA i = 0 hasta N-1:
      // Vecinos con condición de Neumann en i=0 e i=N-1
      i_left  = (i == 0)   ? 1     : i - 1
      i_right = (i == N-1) ? N - 2 : i + 1
      
      laplaciano = u_curr[i_right][j] + u_curr[i_left][j]
                 + u_curr[i][j+1]     + u_curr[i][j-1]
                 - 4 * u_curr[i][j]
      
      u_next[i][j] = 2*u_curr[i][j] - u_prev[i][j] + r2 * laplaciano

  // PASO 2: Aplicar condición de Dirichlet en j = 0 (borde fijo)
  PARA i = 0 hasta N-1:
    u_next[i][0] = 0.0

  // PASO 3: Aplicar condición de Dirichlet en j = N-1 (borde forzado)
  PARA i = 0 hasta N-1:
    u_next[i][N-1] = A * sin(omega * t_k)

  // PASO 4: Avanzar en el tiempo
  u_prev = u_curr
  u_curr = u_next

FIN BUCLE
```

El **PASO 1** (el doble bucle sobre todos los nodos) concentra prácticamente la totalidad del coste computacional. Para $N = 1000$ nodos por lado y $10^5$ pasos temporales, supone del orden de $10^{11}$ operaciones en punto flotante. Es precisamente en este bucle doble donde se aplican las directivas `#pragma omp parallel for` de OpenMP para distribuir las iteraciones entre los hilos disponibles, logrando una reducción del tiempo de ejecución proporcional al número de núcleos físicos utilizados.
```
