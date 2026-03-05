# Simulación de los Modos Normales de Vibración de una Sábana Elástica mediante Diferencias Finitas

## Tabla de contenidos

1. [Introducción y motivación](#1-introducción-y-motivación)
2. [Descripción física del sistema](#2-descripción-física-del-sistema)
3. [Modelo físico completo: la ecuación gobernante](#3-modelo-físico-completo-la-ecuación-gobernante)
4. [Condiciones de contorno](#4-condiciones-de-contorno)
5. [Perfil de equilibrio estático](#5-perfil-de-equilibrio-estático)
6. [Soluciones analíticas: modos normales de vibración](#6-soluciones-analíticas-modos-normales-de-vibración)
7. [Respuesta en frecuencia y resonancia amortiguada](#7-respuesta-en-frecuencia-y-resonancia-amortiguada)
8. [Discretización del dominio espacial](#8-discretización-del-dominio-espacial)
9. [Discretización temporal: esquema en diferencias finitas](#9-discretización-temporal-esquema-en-diferencias-finitas)
10. [Discretización de las condiciones de contorno](#10-discretización-de-las-condiciones-de-contorno)
11. [Condición de estabilidad de Courant-Friedrichs-Lewy (CFL)](#11-condición-de-estabilidad-de-courant-friedrichs-lewy-cfl)
12. [Resumen del algoritmo de integración temporal](#12-resumen-del-algoritmo-de-integración-temporal)
13. [Efectos físicos no modelados](#13-efectos-físicos-no-modelados)

---

## 1. Introducción y motivación

Una sábana extendida y sometida a tensión es un sistema físico continuo cuya dinámica está gobernada por una ecuación en derivadas parciales (EDP) hiperbólica: la **ecuación de onda bidimensional**. Sin embargo, una sábana real no se comporta como una membrana ideal: el material textil es **viscoelástico** (combina propiedades elásticas y viscosas internas), el aire que la rodea introduce **amortiguamiento proporcional a la velocidad**, y la gravedad actúa como una **carga transversal distribuida** que deforma el perfil de reposo.

Cuando uno de sus bordes es forzado a oscilar a una frecuencia determinada, la sábana responde con una superposición de sus modos propios de vibración. Si la frecuencia de forzamiento coincide con alguna de las frecuencias naturales del sistema se produce **resonancia**: la amplitud de oscilación alcanza un máximo **finito**, determinado por el amortiguamiento total del sistema. Este comportamiento contrasta con el de la membrana ideal, donde la amplitud en resonancia crecería indefinidamente, y es exactamente lo que se observa en la realidad.

Desde el punto de vista computacional, la integración numérica de la ecuación de onda 2D sobre una malla fina constituye un problema de alto coste computacional: en cada paso de tiempo es necesario actualizar el valor de la función en los $N^2$ puntos de la malla, y el coste total escala como $\mathcal{O}(N^3)$ para un tiempo de simulación fijo. Este carácter intensivo en cómputo lo convierte en un caso de uso ideal para la **paralelización con OpenMP**, cuyo impacto en el rendimiento es el objeto principal de este trabajo.

---

## 2. Descripción física del sistema

Se considera una sábana delgada de forma cuadrada con lado de longitud $L$, dispuesta **horizontalmente** en el plano $XY$ y sometida a una tensión superficial uniforme $T$ (en N/m). El desplazamiento se produce en la dirección **vertical** $Z$, que es también la dirección de la gravedad.

### 2.1 Sistema de referencia

- El plano $XY$ define el plano de reposo de la sábana.
- El dominio espacial es el cuadrado $\Omega = [0, L] \times [0, L]$.
- El desplazamiento transversal en el punto $(x, y)$ y en el instante $t$ se denota $u(x, y, t)$ [m], tomado positivo hacia arriba.
- La gravedad apunta en la dirección $-Z$, con módulo $g = 9.81$ m/s².

### 2.2 Hipótesis del modelo

- Los desplazamientos transversales $u$ son **pequeños** en comparación con $L$, de modo que la tensión superficial puede considerarse constante, uniforme e isótropa durante el movimiento (hipótesis de membrana linealizada).
- La densidad superficial de masa $\rho$ [kg/m²] es uniforme.
- El material es **viscoelástico lineal**, modelado mediante el modelo de Kelvin-Voigt: la respuesta mecánica del tejido combina un componente elástico (proporcional a la deformación) y un componente viscoso interno (proporcional a la velocidad de deformación), actuando en paralelo.
- El rozamiento con el aire es **viscoso lineal**: la fuerza de amortiguamiento por unidad de área es proporcional a la velocidad local $\partial u / \partial t$ y se opone al movimiento.
- La tensión superficial es la misma en todas las direcciones del plano (material isótropo).

### 2.3 Configuración de los bordes

| Borde | Ecuación del contorno | Condición física |
|---|---|---|
| Borde inferior | $y = 0$ | **Fijo** — nodo permanente, anclado a una estructura rígida |
| Borde superior | $y = L$ | **Forzado** — movimiento armónico prescrito en $Z$ |
| Borde izquierdo | $x = 0$ | **Libre** — sin fuerza normal externa |
| Borde derecho | $x = L$ | **Libre** — sin fuerza normal externa |

---

## 3. Modelo físico completo: la ecuación gobernante

### 3.1 Fuerzas que actúan sobre un elemento diferencial

Consideremos un elemento diferencial de sábana de dimensiones $dx \times dy$, masa $dm = \rho \, dx \, dy$, y desplazamiento vertical $u(x,y,t)$. Sobre este elemento actúan cuatro tipos de fuerzas en la dirección $Z$:

**a) Fuerza elástica de membrana (tensión superficial)**

La tensión superficial $T$ actúa tangencialmente a lo largo de los cuatro bordes del elemento diferencial. Para desplazamientos pequeños, la componente vertical neta de estas fuerzas de tensión es proporcional a la curvatura local de la superficie:

$$dF_\text{elástica} = T \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) dx \, dy = T \, \nabla^2 u \, dx \, dy$$

El operador laplaciano $\nabla^2 u$ mide la curvatura media de la superficie: si la sábana está curvada hacia arriba en un punto, la fuerza de tensión actúa hacia abajo restaurando el equilibrio, y viceversa.

**b) Fuerza viscosa interna del material: modelo de Kelvin-Voigt**

En el modelo de Kelvin-Voigt, el material textil se idealiza como un resorte elástico y un amortiguador viscoso dispuestos **en paralelo**. El resorte reproduce la respuesta elástica (ya capturada por el término de tensión superficial). El amortiguador interno genera una fuerza proporcional a la velocidad de deformación: cuanto más rápido cambia la curvatura de la sábana, mayor es la fuerza viscosa interna que se opone a ese cambio. Matemáticamente, la fuerza viscosa interna por unidad de área en la dirección $Z$ es:

$$dF_\text{visc. interna} = \eta \, \nabla^2 \dot{u} \, dx \, dy = \eta \left( \frac{\partial^2 \dot{u}}{\partial x^2} + \frac{\partial^2 \dot{u}}{\partial y^2} \right) dx \, dy$$

donde $\dot{u} = \partial u / \partial t$ es la velocidad transversal y $\eta$ [N·s/m] es el **coeficiente de viscosidad superficial interna** del material. La consecuencia física más importante de este término es que **amortigua preferentemente los modos de alta frecuencia**: los modos con muchos nodos tienen mayor curvatura (mayor laplaciano) y por tanto experimentan mayor fuerza viscosa interna, disipándose mucho más rápido que los modos bajos.

**c) Amortiguamiento viscoso externo (rozamiento con el aire)**

El aire en contacto con ambas caras de la sábana ejerce sobre ella una fuerza de rozamiento proporcional a la velocidad local y opuesta al movimiento:

$$dF_\text{aire} = -b \, \dot{u} \, dx \, dy$$

donde $b$ [kg/(m²·s)] es el **coeficiente de amortiguamiento viscoso externo**. A diferencia del amortiguamiento viscoelástico interno, este término amortigua todos los modos de forma aproximadamente **uniforme**, independientemente de su frecuencia espacial.

**d) Gravedad**

Al estar la sábana horizontal y las oscilaciones en la dirección vertical, la gravedad actúa como una carga transversal uniforme y constante sobre toda la membrana:

$$dF_\text{gravedad} = -\rho \, g \, dx \, dy$$

En esta configuración la gravedad no modifica la tensión superficial en el plano (como sí haría en una membrana vertical), sino que simplemente actúa como término fuente: su efecto neto es **combrar la sábana hacia abajo** en su posición de reposo.

### 3.2 Ecuación de movimiento completa

Aplicando la segunda ley de Newton al elemento diferencial y sumando todas las contribuciones:

$$\rho \, \frac{\partial^2 u}{\partial t^2} = T \, \nabla^2 u + \eta \, \nabla^2 \dot{u} - b \, \dot{u} - \rho \, g$$

Definiendo los parámetros reducidos del modelo:

| Parámetro | Definición | Unidades | Significado físico |
|---|---|---|---|
| $c$ | $\sqrt{T/\rho}$ | m/s | Velocidad de propagación de ondas elásticas |
| $\alpha$ | $\eta / \rho$ | m²/s | Difusividad viscoelástica del material |
| $\beta$ | $b / \rho$ | 1/s | Tasa de amortiguamiento viscoso externo |

La **ecuación gobernante** en su forma canónica es:

$$\boxed{\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u + \alpha \, \nabla^2 \dot{u} - \beta \, \dot{u} - g}$$

### 3.3 Interpretación física de cada término

- $c^2 \nabla^2 u$: fuerza restauradora elástica. Término dominante; gobierna la propagación de ondas transversales con velocidad $c$.
- $\alpha \, \nabla^2 \dot{u}$: disipación viscoelástica interna. Actúa como una difusión de velocidades y amortigua selectivamente los modos de alto número de onda.
- $-\beta \, \dot{u}$: amortiguamiento por rozamiento con el aire. Frena el movimiento uniformemente en toda la sábana.
- $-g$: carga gravitatoria uniforme. Desplaza el perfil de equilibrio hacia abajo sin alterar la dinámica de las perturbaciones.

### 3.4 Descomposición en equilibrio y perturbación

Es conveniente descomponer el desplazamiento en:

$$u(x, y, t) = u_\text{eq}(x, y) + u'(x, y, t)$$

El perfil de equilibrio estático satisface la ecuación de Poisson $c^2 \nabla^2 u_\text{eq} = g$, y la perturbación dinámica satisface la ecuación sin el término gravitatorio. En la simulación numérica se integra directamente la ecuación completa, capturando ambos efectos de forma automática.

---

## 4. Condiciones de contorno

### 4.1 Borde fijo: $y = 0$ (Dirichlet homogéneo)

$$u(x, 0, t) = 0 \quad \forall x \in [0, L], \, \forall t \geq 0$$

### 4.2 Borde forzado: $y = L$ (Dirichlet no homogéneo)

$$u(x, L, t) = A \sin(\omega t) \quad \forall x \in [0, L], \, \forall t \geq 0$$

El forzamiento es uniforme en $x$, por lo que inyecta energía principalmente en los modos con $m=0$ (perfil constante en $x$). Cuando su frecuencia coincide con alguna $\omega_{mn}$, el sistema entra en resonancia.

### 4.3 Bordes libres: $x = 0$ y $x = L$ (Neumann homogéneo)

$$\frac{\partial u}{\partial x}(0, y, t) = 0, \quad \frac{\partial u}{\partial x}(L, y, t) = 0 \quad \forall y \in [0, L], \, \forall t \geq 0$$

La sábana puede adoptar cualquier altura en los bordes laterales, pero su tangente es siempre horizontal en la dirección $x$ en esos puntos. Las ondas se reflejan perfectamente en estos bordes sin pérdida de energía.

### 4.4 Condiciones iniciales

La sábana parte del reposo en su posición no deformada:

$$u(x, y, 0) = 0, \qquad \frac{\partial u}{\partial t}(x, y, 0) = 0$$

---

## 5. Perfil de equilibrio estático

El perfil de equilibrio $u_\text{eq}(y)$ satisface $c^2 d^2u_\text{eq}/dy^2 = g$ con $u_\text{eq}(0)=0$ y $u_\text{eq}(L)=0$. Integrando dos veces:

$$\boxed{u_\text{eq}(y) = \frac{g}{2c^2} \, y(y - L)}$$

Esta es una parábola cóncava hacia abajo con descenso máximo en el centro:

$$u_\text{eq,min} = u_\text{eq}(L/2) = -\frac{g L^2}{8c^2}$$

Para una sábana típica de $L = 2$ m, $T = 20$ N/m, $\rho = 0.2$ kg/m² (con $c = 10$ m/s), el descenso central es $u_\text{eq,min} \approx -0.49$ mm: pequeño pero no despreciable para amplitudes de oscilación del orden del milímetro.

---

## 6. Soluciones analíticas: modos normales de vibración

### 6.1 Separación de variables

Para el sistema sin forzamiento ni amortiguamiento ($\alpha = \beta = 0$, $g = 0$), se buscan soluciones de la forma $u' = X(x)Y(y)T(t)$. Las condiciones de contorno seleccionan:

**En $x$ (condiciones de Neumann):**
$$X_m(x) = \cos\left(\frac{m\pi x}{L}\right), \quad m = 0, 1, 2, \ldots$$

**En $y$ (condiciones de Dirichlet):**
$$Y_n(y) = \sin\left(\frac{n\pi y}{L}\right), \quad n = 1, 2, 3, \ldots$$

### 6.2 Formas propias y frecuencias naturales

Las **formas propias** son:

$$\phi_{mn}(x,y) = \cos\left(\frac{m\pi x}{L}\right) \sin\left(\frac{n\pi y}{L}\right)$$

y las **frecuencias angulares naturales**:

$$\boxed{\omega_{mn} = \frac{c\pi}{L}\sqrt{m^2 + n^2}}$$

Los primeros modos ordenados por frecuencia creciente:

| Modo $(m, n)$ | $\sqrt{m^2+n^2}$ | $f_{mn}/f_{01}$ | Descripción |
|:---:|:---:|:---:|---|
| $(0, 1)$ | $1$ | $1.000$ | Fundamental: uniforme en $x$, una semionda en $y$ |
| $(1, 1)$ | $\sqrt{2} \approx 1.414$ | $1.414$ | Una semionda en $x$, una en $y$ |
| $(0, 2)$ | $2$ | $2.000$ | Uniforme en $x$, dos semiondes en $y$ |
| $(1, 2)$ | $\sqrt{5} \approx 2.236$ | $2.236$ | Una en $x$, dos en $y$ |
| $(2, 1)$ | $\sqrt{5} \approx 2.236$ | $2.236$ | Dos en $x$, una en $y$ |
| $(2, 2)$ | $2\sqrt{2} \approx 2.828$ | $2.828$ | Dos en $x$, dos en $y$ |
| $(0, 3)$ | $3$ | $3.000$ | Uniforme en $x$, tres semiondes en $y$ |

---

## 7. Respuesta en frecuencia y resonancia amortiguada

### 7.1 Ecuación modal con amortiguamiento

Al proyectar la ecuación gobernante sobre cada modo propio $\phi_{mn}$ se obtiene la ecuación del oscilador armónico amortiguado:

$$\ddot{q}_{mn} + 2\zeta_{mn}\omega_{mn}\dot{q}_{mn} + \omega_{mn}^2 q_{mn} = F_{mn}(t)$$

El **factor de amortiguamiento modal** recibe contribuciones de ambos mecanismos de disipación:

$$\boxed{\zeta_{mn} = \underbrace{\frac{\beta}{2\omega_{mn}}}_{\text{amort. externo (uniforme)}} + \underbrace{\frac{\alpha \, k_{mn}^2}{2\omega_{mn}}}_{\text{amort. viscoelástico (crece con el modo)}}}$$

con $k_{mn}^2 = (m^2 + n^2)\pi^2/L^2$.

Este resultado formaliza la diferencia entre los dos mecanismos: el amortiguamiento externo $\beta$ contribuye de forma que decrece con la frecuencia del modo, mientras que el viscoelástico $\alpha k_{mn}^2$ crece con $\sqrt{m^2+n^2}$, disipando mucho más rápido los modos altos. La combinación produce un amortiguamiento no monótono en función del número de modo, que es físicamente realista.

### 7.2 Amplitud en resonancia

Para forzamiento armónico $F_{mn}(t) = F_0\sin(\omega t)$, la amplitud en régimen estacionario es:

$$|q_{mn}|_\text{estac} = \frac{F_0/\omega_{mn}^2}{\sqrt{(1-\Omega^2)^2 + (2\zeta_{mn}\Omega)^2}}, \quad \Omega = \frac{\omega}{\omega_{mn}}$$

En resonancia exacta ($\Omega = 1$):

$$|q_{mn}|_\text{res} = \frac{F_0/\omega_{mn}^2}{2\zeta_{mn}}$$

La amplitud es **finita** siempre que $\zeta_{mn} > 0$, a diferencia del sistema ideal donde crecería linealmente con el tiempo. El factor de calidad $Q = 1/(2\zeta_{mn})$ mide la agudeza del pico: para $\zeta_{mn} = 0.01$ (amortiguamiento del 1%), $Q = 50$ y la amplitud en resonancia es 50 veces la amplitud estática.

---

## 8. Discretización del dominio espacial

### 8.1 Malla uniforme

El dominio $\Omega = [0,L]\times[0,L]$ se discretiza con una malla uniforme de $N \times N$ nodos y paso espacial:

$$h = \frac{L}{N-1}, \quad x_i = ih, \quad y_j = jh, \quad i,j = 0,\ldots,N-1$$

Se denota $u_{i,j}^k \approx u(x_i, y_j, t_k)$ con $t_k = k\Delta t$.

### 8.2 Laplaciano discreto: stencil de 5 puntos

Por desarrollo de Taylor de orden 4 se obtiene la aproximación de segundo orden:

$$\left(\nabla^2 u\right)_{i,j}^k \approx \frac{u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k}{h^2} + \mathcal{O}(h^2)$$

El error de truncación dominante es $-(h^2/12)(\partial^4u/\partial x^4 + \partial^4u/\partial y^4)$, confirmando el orden 2 espacial.

### 8.3 Discretización del término viscoelástico

Se define el **numerador del laplaciano discreto**:

$$\mathcal{L}_{i,j}^k = u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k$$

La velocidad se aproxima con diferencias hacia atrás para mantener la explicitud del esquema:

$$\dot{u}_{i,j}^k \approx \frac{u_{i,j}^k - u_{i,j}^{k-1}}{\Delta t} + \mathcal{O}(\Delta t)$$

El laplaciano de la velocidad queda entonces:

$$\left(\nabla^2\dot{u}\right)_{i,j}^k \approx \frac{\mathcal{L}_{i,j}^k - \mathcal{L}_{i,j}^{k-1}}{h^2 \Delta t}$$

---

## 9. Discretización temporal: esquema en diferencias finitas

### 9.1 Ensamblaje de la ecuación discreta completa

Aplicando diferencias centradas de segundo orden para la aceleración y sustituyendo todas las aproximaciones en la ecuación gobernante:

$$\frac{u_{i,j}^{k+1} - 2u_{i,j}^k + u_{i,j}^{k-1}}{\Delta t^2} = \frac{c^2}{h^2}\mathcal{L}_{i,j}^k + \frac{\alpha}{h^2\Delta t}\left(\mathcal{L}_{i,j}^k - \mathcal{L}_{i,j}^{k-1}\right) - \beta\frac{u_{i,j}^k - u_{i,j}^{k-1}}{\Delta t} - g$$

Despejando $u_{i,j}^{k+1}$ y definiendo los coeficientes adimensionales $r^2 = c^2\Delta t^2/h^2$ y $\gamma = \alpha\Delta t/h^2$:

$$\boxed{u_{i,j}^{k+1} = (2-\beta\Delta t)\,u_{i,j}^k + (\beta\Delta t - 1)\,u_{i,j}^{k-1} + (r^2+\gamma)\,\mathcal{L}_{i,j}^k - \gamma\,\mathcal{L}_{i,j}^{k-1} - g\Delta t^2}$$

Para $\alpha=0$, $\beta=0$, $g=0$ se recupera exactamente el esquema clásico de la membrana ideal: $u_{i,j}^{k+1} = 2u_{i,j}^k - u_{i,j}^{k-1} + r^2\mathcal{L}_{i,j}^k$.

### 9.2 Propiedades del esquema

| Propiedad | Valor |
|---|---|
| Tipo | Explícito (leapfrog modificado) |
| Orden temporal | $\mathcal{O}(\Delta t^2)$ elástico/gravitatorio; $\mathcal{O}(\Delta t)$ disipativos |
| Orden espacial | $\mathcal{O}(h^2)$ |
| Coste por paso | $\mathcal{O}(N^2)$ operaciones |
| Sistema lineal a resolver | No |
| Arrays necesarios | 3 niveles de $u$ + laplaciano del nivel anterior $\mathcal{L}^{k-1}$ |

### 9.3 Inicialización del primer paso

Con condición inicial de velocidad nula, la simetría $u^{-1} = u^1$ permite deducir:

$$u_{i,j}^1 = -\frac{g\Delta t^2}{2}$$

para los nodos interiores (la sábana comienza a caer por la gravedad desde el primer instante).

---

## 10. Discretización de las condiciones de contorno

### 10.1 Borde fijo ($j = 0$)

$$u_{i,0}^k = 0 \quad \forall i, k$$

### 10.2 Borde forzado ($j = N-1$)

$$u_{i,N-1}^k = A\sin(\omega\,k\,\Delta t) \quad \forall i$$

Esta condición se impone después de la actualización general, sobreescribiendo los valores calculados en esa fila.

### 10.3 Bordes libres: nodos fantasma

La condición $\partial u/\partial x = 0$ en $x=0$ y $x=L$ se implementa con nodos virtuales fuera del dominio. Imponiendo diferencias centradas de orden 2:

$$u_{-1,j}^k = u_{1,j}^k, \quad u_{N,j}^k = u_{N-2,j}^k$$

El laplaciano discreto en los nodos del borde queda:

$$\mathcal{L}_{0,j}^k = 2u_{1,j}^k + u_{0,j+1}^k + u_{0,j-1}^k - 4u_{0,j}^k$$

$$\mathcal{L}_{N-1,j}^k = 2u_{N-2,j}^k + u_{N-1,j+1}^k + u_{N-1,j-1}^k - 4u_{N-1,j}^k$$

La misma corrección se aplica a $\mathcal{L}^{k-1}$ de forma consistente.

---

## 11. Condición de estabilidad de Courant-Friedrichs-Lewy (CFL)

### 11.1 Condición hiperbólica (onda elástica)

El análisis de estabilidad de Von Neumann para el término elástico, considerando el caso más restrictivo (propagación diagonal), conduce a:

$$r = \frac{c\Delta t}{h} \leq \frac{1}{\sqrt{2}} \approx 0.707$$

### 11.2 Condición parabólica (viscoelasticidad)

El término $\alpha\nabla^2\dot{u}$ introduce comportamiento difusivo con la condición adicional:

$$\gamma = \frac{\alpha\Delta t}{h^2} \leq \frac{1}{4}$$

### 11.3 Condición de diseño

Ambas deben satisfacerse simultáneamente. En la práctica, para materiales textiles reales, la condición hiperbólica es siempre la más restrictiva. Se adopta $r = 0.5$ como valor de diseño:

$$\boxed{\Delta t = \frac{0.5\,h}{c}}$$

Con este $\Delta t$, la condición parabólica exige $\alpha \leq ch/2$, lo que debe verificarse explícitamente al fijar los parámetros del material.

### 11.4 Implicaciones en el coste computacional

La condición CFL impone $\Delta t \propto 1/N$, de modo que al doblar $N$: el número de nodos se cuadruplica y el número de pasos se dobla. El coste total escala como $\mathcal{O}(N^3)$, lo que justifica plenamente la paralelización con OpenMP.

---

## 12. Resumen del algoritmo de integración temporal

```
PARÁMETROS FÍSICOS DE ENTRADA:
  L       [m]         Longitud del lado de la sábana
  rho     [kg/m²]     Densidad superficial de masa
  T_tens  [N/m]       Tensión superficial
  eta     [N·s/m]     Coeficiente de viscosidad interna (Kelvin-Voigt)
  b       [kg/(m²s)]  Coeficiente de amortiguamiento viscoso externo
  A       [m]         Amplitud del forzamiento
  omega   [rad/s]     Frecuencia angular del forzamiento
  T_sim   [s]         Tiempo total de simulación
  N       [-]         Número de nodos por lado

PREPROCESO Y VERIFICACIÓN:
  c     = sqrt(T_tens / rho)
  alpha = eta / rho
  beta  = b / rho
  h     = L / (N - 1)
  dt    = 0.5 * h / c               // CFL hiperbólica (r = 0.5)
  r2    = (c * dt / h)^2            // = 0.25
  gam   = alpha * dt / h^2          // coeficiente viscoelástico
  ASSERT gam <= 0.25                 // verificar CFL parabólica
  Nstep = ceil(T_sim / dt)

INICIALIZACIÓN:
  u_prev[i][j] = 0     para todo (i,j)
  u_curr[i][j] = 0     para todo (i,j)
  L_prev[i][j] = 0     para todo (i,j)

BUCLE TEMPORAL (k = 1, 2, ..., Nstep):

  t_k = k * dt

  // PASO 1: Calcular laplaciano L_curr con condiciones de Neumann
  PARA j = 1 hasta N-2:
    PARA i = 0 hasta N-1:
      i_l = (i == 0)   ? 1   : i-1
      i_r = (i == N-1) ? N-2 : i+1
      L_curr[i][j] = u_curr[i_r][j] + u_curr[i_l][j]
                   + u_curr[i][j+1] + u_curr[i][j-1]
                   - 4.0 * u_curr[i][j]

  // PASO 2: Actualizar todos los nodos interiores
  PARA j = 1 hasta N-2:
    PARA i = 0 hasta N-1:
      u_next[i][j] = (2.0 - beta*dt) * u_curr[i][j]
                   + (beta*dt - 1.0) * u_prev[i][j]
                   + (r2 + gam)      * L_curr[i][j]
                   -  gam            * L_prev[i][j]
                   -  g * dt * dt

  // PASO 3: Imponer CC en borde fijo (j = 0)
  PARA i = 0 hasta N-1:
    u_next[i][0] = 0.0

  // PASO 4: Imponer CC en borde forzado (j = N-1)
  PARA i = 0 hasta N-1:
    u_next[i][N-1] = A * sin(omega * t_k)

  // PASO 5: Rotar arrays temporales
  L_prev = L_curr
  u_prev = u_curr
  u_curr = u_next

FIN BUCLE
```

El núcleo computacional son los pasos 1 y 2: dos bucles dobles sobre los $N \times (N-2)$ nodos. Las iteraciones del bucle externo (sobre $j$) son completamente independientes entre sí, ya que la actualización de cada nodo utiliza únicamente valores del instante anterior. Esta independencia de datos es la que hace posible la paralelización directa con `#pragma omp parallel for` sin ninguna condición de carrera.

---

## 13. Efectos físicos no modelados

| Efecto | Descripción | Impacto estimado |
|---|---|---|
| **Rigidez a la flexión** | La tela resiste curvarse (término $D\nabla^4 u$). Relevante para telas gruesas. | Modifica frecuencias propias de modos altos |
| **Anisotropía del material** | Urdimbre y trama con propiedades distintas: $T_x \neq T_y$. | Rompe la simetría de los modos; separa modos degenerados |
| **No linealidades geométricas** | Para $u \sim L$, la tensión ya no es constante. | Acopla modos; genera frecuencias de combinación |
| **Amortiguamiento aerodinámico no lineal** | Rozamiento $\propto \dot{u}^2$ a velocidades altas. | Solo relevante para oscilaciones muy rápidas |
| **Masa añadida del aire** | El aire que se mueve con la sábana aumenta la inercia efectiva. | Reduce ligeramente todas las frecuencias propias |
| **No uniformidad real de la tensión** | Depende de cómo se fija la sábana en los bordes. | Difícil de cuantificar sin medición experimental |

Para los objetivos de este trabajo, el modelo adoptado representa el punto óptimo entre **fidelidad física** y **tractabilidad numérica y computacional**.
