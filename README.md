# 🦠 Modelo SIR — Simulador Epidemiológico

![HTML5](https://img.shields.io/badge/HTML5-E34F26?style=flat-square&logo=html5&logoColor=white)
![CSS3](https://img.shields.io/badge/CSS3-1572B6?style=flat-square&logo=css3&logoColor=white)
![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?style=flat-square&logo=javascript&logoColor=black)
![Plotly](https://img.shields.io/badge/Plotly.js-3F4F75?style=flat-square&logo=plotly&logoColor=white)
![Ecuaciones Diferenciales](https://img.shields.io/badge/Ecuaciones%20Diferenciales-2025--2026-58a6ff?style=flat-square)

Simulador web del modelo SIR (Susceptibles–Infectados–Recuperados) para la visualización interactiva de brotes epidémicos en entornos universitarios. Proyecto final de la asignatura **Ecuaciones Diferenciales**.

---

## 📖 Descripción

El modelo SIR es un sistema de ecuaciones diferenciales ordinarias que describe cómo se propaga una enfermedad infecciosa en una población cerrada. Este simulador permite explorar visualmente el comportamiento del modelo ajustando parámetros en tiempo real, sin necesidad de instalar nada — solo un navegador.

El objetivo es **didáctico**: mostrar cómo las ecuaciones diferenciales modelan fenómenos reales, comparar métodos numéricos de resolución (Euler, Heun, RK4), y entender conceptos como el número reproductivo básico R₀ o el efecto de la vacunación.

Desarrollado como proyecto final para la asignatura **Ecuaciones Diferenciales**, semestre 2025–2026.

---

## ✨ Funcionalidades

- **🧮 Tres métodos numéricos** — Euler (orden 1), Heun (orden 2) y Runge-Kutta RK4 (orden 4) implementados manualmente y comparados en una misma gráfica con tabla de errores relativos.
- **🎚️ Sliders interactivos** — 7 parámetros ajustables (N, I₀, β, γ, vacunación, duración, dt) con actualización instantánea de la gráfica.
- **📊 Estadísticas clave** — R₀ con código de color dinámico (verde/amarillo/rojo), día del pico, máximo de infectados simultáneos, total de recuperados y porcentaje afectado.
- **🎴 Escenarios preconfigurados** — 4 cards cliqueables: Semana de parciales, Campaña de vacunación, Modalidad virtual y Brote descontrolado.
- **📈 Derivadas en vivo** — Las ecuaciones diferenciales se muestran con los valores numéricos sustituidos y el resultado calculado en tiempo real, coloreado según el signo.
- **🔬 Modo paso a paso** — Avanza día por día sobre la solución mostrando los valores actuales de S, I, R, las derivadas y las 4 pendientes intermedias del método RK4 (k₁, k₂, k₃, k₄).
- **🖱️ Ecuaciones cliqueables** — Al hacer clic en una ecuación (dS/dt, dI/dt, dR/dt), la curva correspondiente en la gráfica se resalta y las otras se atenúan.
- **⏱️ Paso de tiempo ajustable** — Slider dt (0.1–5.0) para modificar la resolución temporal y observar cómo afecta la precisión de los métodos numéricos.
- **📐 Comparación de escenarios** — Las cards permiten saltar entre configuraciones predefinidas sin manipular sliders manualmente.
- **🌙 Tema oscuro** — Interfaz científica con paleta oscura, tipografía monoespaciada para ecuaciones y diseño responsivo.

---

## 🧮 Modelo Matemático

El sistema de ecuaciones diferenciales del modelo SIR con **incidencia estándar** (división por N):

```math
\frac{dS}{dt} = -\beta \frac{S I}{N}
```

```math
\frac{dI}{dt} = \beta \frac{S I}{N} - \gamma I
```

```math
\frac{dR}{dt} = \gamma I
```

| Símbolo | Significado | Unidad |
|---|---|---|
| S(t) | Población susceptible | individuos |
| I(t) | Población infectada activa | individuos |
| R(t) | Población recuperada / inmune | individuos |
| N | Población total (S + I + R) | individuos |
| β | Tasa de contagio (contactos efectivos por día) | contactos/día |
| γ | Tasa de recuperación | 1/días |
| R₀ | Número reproductivo básico = β / γ | adimensional |

La división por **N** normaliza la probabilidad de contacto, haciendo que β sea una propiedad intrínseca de la enfermedad independiente del tamaño poblacional. Si R₀ > 1 la epidemia crece; si R₀ < 1 se extingue.

---

## 🔢 Métodos Numéricos

La aplicación resuelve el sistema con tres métodos, todos implementados manualmente en JavaScript:

| Método | Orden | Descripción | Precisión |
|---|---|---|---|
| Euler | 1 | Aproximación lineal con la derivada actual | Baja (error O(dt)) |
| Heun (Euler mejorado) | 2 | Predictor-corrector con promedio de dos pendientes | Media (error O(dt²)) |
| Runge-Kutta RK4 | 4 | Cuatro evaluaciones de pendiente por paso | Alta (error O(dt⁴)) |

La gráfica de comparación numérica superpone las curvas de I(t) de los tres métodos y calcula el **error relativo máximo** de Euler y Heun respecto a RK4, que se toma como referencia.

---

## 🖥️ Guía de la Interfaz

La aplicación tiene dos pestañas navegables desde la barra superior.

### 📖 Pestaña: Fundamentos

Sección informativa del modelo SIR:
- **Compartimentos** — Cards con los 3 grupos poblacionales (S, I, R) y diagrama de flujo S → I → R.
- **Ecuaciones** — Las EDOs del modelo con nota explicativa sobre la normalización por N y chips de parámetros (β, γ, R₀).
- **Interpretación** — Significado físico de cada derivada.
- **Supuestos** — Lista de hipótesis del modelo (población constante, mezcla homogénea, inmunidad permanente).
- **Aplicaciones** — Usos del modelo en salud pública, vacunación e investigación.
- **Cómo usar** — 5 pasos para operar el simulador.

### 📊 Pestaña: Simulación

Panel de control interactivo y visualización:

| Componente | Descripción |
|---|---|
| **Sliders** | N (1000–50000), I₀ (1–100), β (0.05–1.00), γ (0.01–0.50), % vacunado (0–80%), duración (30–365 días), dt (0.1–5.0). Cada slider muestra su valor actual en un badge. |
| **Estadísticas clave** | R₀ con código de color, día del pico, máximo de infectados, total recuperados, porcentaje afectado. Se actualizan en tiempo real. |
| **Gráfica principal** | Curvas S(t) en azul, I(t) en rojo, R(t) en verde. Línea vertical punteada marcando el pico de infectados. Tooltips al hover con valores exactos. |
| **Escenarios** | 4 cards que reconfiguran los sliders automáticamente y ejecutan la simulación. |
| **Comparación numérica** | Segunda gráfica con I(t) de Euler, Heun y RK4 superpuestas. Tabla con error relativo máximo. |
| **Derivadas en vivo** | Muestra la evaluación numérica de cada derivada con los valores actuales sustituidos y el resultado coloreado (verde si es positivo, rojo si es negativo). |
| **Modo paso a paso** | Toggle switch + navegación día por día. Muestra S, I, R actuales, las 3 derivadas y la tabla completa de pendientes RK4 (k₁–k₄). Un cursor azul punteado se desplaza sobre la gráfica. |
| **Ecuaciones cliqueables** | Al hacer clic en una ecuación del panel superior, la curva correspondiente en la gráfica se resalta y las demás se atenúan. |

---

## 🗂️ Estructura del Proyecto

```
proyecto/
├── index.html        — Estructura HTML5 con ambas pestañas y todo el markup
├── simulation.js     — Lógica: solvers, estado, gráficas, eventos (~840 líneas)
└── styles.css        — Tema oscuro responsivo (~1180 líneas)
```

Sin frameworks, sin bundlers, sin backend. Todo funciona abriendo `index.html` directamente en el navegador (protocolo `file://`).

---

## 🚀 Cómo Usar

1. **Abrir** `index.html` en cualquier navegador moderno (Chrome, Firefox, Edge).
2. **Ajustar** los sliders del panel de control para modificar los parámetros de la simulación.
3. **Explorar** los escenarios preconfigurados haciendo clic en las cards.
4. **Activar** el modo paso a paso para observar la evolución día por día y las pendientes del método RK4.

No requiere instalación, servidor ni conexión a internet (excepto la primera carga de Plotly.js desde CDN).

---

## 🛠️ Tecnologías

- **HTML5** — Estructura semántica y markup de la interfaz
- **CSS3** — Diseño oscuro, grid layout, animaciones y responsividad
- **JavaScript** — Lógica de simulación, solvers numéricos y manipulación del DOM
- **Plotly.js** v2.26.0 — Visualización de gráficas interactivas (cargado desde CDN)

---

## 📚 Referencias

- Kermack, W. O. & McKendrick, A. G. (1927). *A contribution to the mathematical theory of epidemics*. Proceedings of the Royal Society A, 115(772), 700–721.
- [Modelo SIR — Wikipedia](https://es.wikipedia.org/wiki/Modelo_SIR)
- [Runge–Kutta methods — Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
- [Plotly.js Documentation](https://plotly.com/javascript/)

---

## 📄 Licencia

Proyecto académico — Ecuaciones Diferenciales · 2025–2026.
