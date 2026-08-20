# Proyecto 1 - Grupo 4

## Metodos Computacionales en Obras Civiles (MCOC)

Analisis de un marco hiperestatico 2D sometido a asentamientos diferenciales en sus apoyos, resuelto con **OpenSeesPy**.

---

## Estructura del Proyecto

```
Proyecto-1---Grupo-4/
├── marco_2d_asentamientos.py   # Script principal de analisis
├── diagramas_esfuerzos.py      # Generacion de diagramas (DMF, DEC, DFA)
├── DMF.png                     # Diagrama de Momentos Flectores
├── DEC.png                     # Diagrama de Esfuerzos Cortantes
├── DFA.png                     # Diagrama de Fuerzas Axiales
├── .venv/                      # Entorno virtual Python 3.12
└── README.md                   # Este archivo
```

---

## Descripcion del Problema

Se analiza un marco compuesto por **una columna (AB)** y **una viga continua (BD)**, con un punto intermedio C donde se mide el desplazamiento vertical. El marco esta sometido a asentamientos (descensos) en sus apoyos, sin cargas externas adicionales.

### Geometria

| Nodo | Posicion (x, y) | Descripcion |
|------|-----------------|-------------|
| 1 (A) | (0.0, 0.0) | Apoyo empotrado |
| 2 (B) | (0.0, 4.0) | Union columna-viga |
| 3 (C) | (4.0, 4.0) | Punto intermedio de la viga |
| 4 (D) | (8.0, 4.0) | Apoyo articulado |

### Condiciones de Borde

- **Nodo 1 (A):** Empotrado (restringido en X, Y y rotacion)
- **Nodo 4 (D):** Articulado (restringido en X y Y, libre en rotacion)

### Propiedades de la Seccion

- Rigidez a flexion: **EI = 1000 tonf*m²**
- Rigidez axial: **EA = infinito** (A = 1.0e12, E = 1.0)

### Asentamientos Impuestos

- **Nodo 1 (A):** -0.02 m (2 cm hacia abajo)
- **Nodo 4 (D):** -0.04 m (4 cm hacia abajo)

---

## Resultados

### Reacciones

| Nodo | Fx (tonf) | Fy (tonf) | Mz (tonf*m) |
|------|-----------|-----------|-------------|
| 1 (A) | +0.256 | +0.085 | +0.341 |
| 4 (D) | -0.256 | -0.085 | 0.000 |

### Momentos Flectores

| Elemento | M_I (tonf*m) | M_J (tonf*m) |
|----------|-------------|-------------|
| Columna A-B | +0.341 | -0.682 |
| Viga B-D | -0.682 | -0.341 |

### Desplazamiento Vertical en Punto C

**dVc = -27.273 mm**

---

## Como Ejecutar

### Requisitos

- Python 3.12 (instalado via winget)
- OpenSeesPy
- Matplotlib

### Pasos

1. Abrir la terminal en VS Code (Ctrl+`)

2. Activar el entorno virtual:
```bash
.venv\Scripts\activate
```

3. Ejecutar el analisis:
```bash
python marco_2d_asentamientos.py
```

4. Generar los diagramas:
```bash
python diagramas_esfuerzos.py
```

---

## Convencion de Signos

- **M_I:** Momento en el extremo I (extremo inicial del elemento), signo negado de OpenSees
- **M_J:** Momento en el extremo J (extremo final del elemento), signo original de OpenSees
- **Momento positivo:** Traccion interna segun convencion de ingenieria estructural
- **Cortante:** Segun sentido local del elemento
- **Axial:** Tension positiva

---

## Tecnologias Utilizadas

- **OpenSeesPy 3.8.0** - Analisis estructural
- **Matplotlib 3.11** - Generacion de graficos
- **Python 3.12** - Lenguaje de programacion
