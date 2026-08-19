"""
Diagramas de Esfuerzos - Marco 2D con Asentamientos
=====================================================
Genera los diagramas de:
1. Momentos Flectores (DMF)
2. Esfuerzos Cortantes (DEC)
3. Fuerzas Axiales (DFA)

Unidades: tonf, m
"""

import matplotlib
matplotlib.use('Agg')
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ============================================================================
# 1. MODELO OPENSEES (mismo que marco_2d_asentamientos.py)
# ============================================================================

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Nodos
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 4.0)
ops.node(3, 4.0, 4.0)
ops.node(4, 8.0, 4.0)

# Apoyos
ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 0)

# Propiedades
A = 1.0e12
E = 1.0
I = 1000.0

# Transformacion y elementos
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, A, E, I, 1)
ops.element('elasticBeamColumn', 2, 2, 3, A, E, I, 1)
ops.element('elasticBeamColumn', 3, 3, 4, A, E, I, 1)

# Asentamientos
ops.constraints('Transformation')
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.sp(1, 2, -0.02)
ops.sp(4, 2, -0.04)

# Analisis
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-10, 10)
ops.algorithm('Linear')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')
ops.analyze(1)
ops.reactions()

# ============================================================================
# 2. EXTRACCION DE FUERZAS INTERNAS
# ============================================================================

elementos = {
    1: {'nodos': (1, 2), 'nombre': 'Columna A-B', 'tipo': 'columna'},
    2: {'nodos': (2, 3), 'nombre': 'Viga B-C',    'tipo': 'viga'},
    3: {'nodos': (3, 4), 'nombre': 'Viga C-D',    'tipo': 'viga'},
}

coord_nodos = {
    1: np.array([0.0, 0.0]),
    2: np.array([0.0, 4.0]),
    3: np.array([4.0, 4.0]),
    4: np.array([8.0, 4.0]),
}

fuerzas = {}
for tag, info in elementos.items():
    f = ops.eleForce(tag)
    fuerzas[tag] = {
        'Fx_I': f[0], 'Fy_I': f[1], 'Mz_I': -f[2],
        'Fx_J': f[3], 'Fy_J': f[4], 'Mz_J': -f[5],
    }

print("Fuerzas internas extraidas:")
for tag, f in fuerzas.items():
    print(f"  Elem {tag}: Fx={f['Fx_I']:.4f},{f['Fx_J']:.4f} | Fy={f['Fy_I']:.4f},{f['Fy_J']:.4f} | M={f['Mz_I']:.4f},{f['Mz_J']:.4f}")

# ============================================================================
# 3. FUNCIONES DE GRAFICACION
# ============================================================================

def obtener_puntos_elemento(tag):
    ni, nj = elementos[tag]['nodos']
    return coord_nodos[ni], coord_nodos[nj]

def rotar_a_local(pI, pJ, valor):
    dx = pJ[0] - pI[0]
    dy = pJ[1] - pI[1]
    L = np.sqrt(dx**2 + dy**2)
    nx = dx / L
    ny = dy / L
    return np.array([-ny, nx]) * valor

def dibujar_elemento_base(ax, pI, pJ, color='#2c3e50', lw=1.5):
    ax.plot([pI[0], pJ[0]], [pI[1], pJ[1]], color=color, linewidth=lw, zorder=2)

def dibujar_diag_funcional(ax, pI, pJ, val_I, val_J, n_pts=50, scale=1.0,
                            fill_color='#3498db', alpha=0.3, line_color='#2980b9',
                            label_max=True, unit='tonf*m'):
    dx = pJ[0] - pI[0]
    dy = pJ[1] - pI[1]
    L = np.sqrt(dx**2 + dy**2)
    tx = np.array([dx, dy]) / L
    nx = np.array([-tx[1], tx[0]])

    s = np.linspace(0, 1, n_pts)
    vals = val_I + (val_J - val_I) * s
    pts_elem = np.array([pI + t * (pJ - pI) for t in s])
    pts_diag = pts_elem + np.outer(vals * scale, nx)

    ax.fill_between(range(n_pts), 0, vals * scale, alpha=0)

    poly_x = list(pts_diag[:, 0]) + [pts_elem[-1, 0]]
    poly_y = list(pts_diag[:, 1]) + [pts_elem[-1, 1]]
    poly_x = [pts_elem[0, 0]] + poly_x
    poly_y = [pts_elem[0, 1]] + poly_y

    ax.fill(poly_x, poly_y, color=fill_color, alpha=alpha, zorder=1)
    ax.plot(pts_diag[:, 0], pts_diag[:, 1], color=line_color, linewidth=2, zorder=3)

    max_idx = np.argmax(np.abs(vals))
    max_val = vals[max_idx]
    if label_max and abs(max_val) > 1e-6:
        pt_label = pts_diag[max_idx]
        offset = nx * 0.3 * np.sign(max_val)
        if abs(offset[0]) < 0.01 and abs(offset[1]) < 0.01:
            offset = nx * 0.4
        ax.annotate(f'{max_val:.3f}',
                     xy=(pt_label[0], pt_label[1]),
                     xytext=(pt_label[0] + offset[0], pt_label[1] + offset[1]),
                     fontsize=8, fontweight='bold', ha='center', va='center',
                     color=line_color,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=line_color, alpha=0.8),
                     zorder=5)

def dibujar_nodos(ax, escala_nodo=0.15):
    for tag, pos in coord_nodos.items():
        ax.plot(pos[0], pos[1], 'o', color='#2c3e50', markersize=6, zorder=4)
        offset = np.array([0.0, -0.35])
        if tag == 1:
            offset = np.array([-0.25, -0.35])
        elif tag == 4:
            offset = np.array([0.25, -0.35])
        ax.annotate(f'N{tag}', xy=(pos[0], pos[1]),
                     xytext=(pos[0] + offset[0], pos[1] + offset[1]),
                     fontsize=9, fontweight='bold', ha='center', va='top', color='#2c3e50')

def dibujar_apoyos(ax):
    sz = 0.3

    p1 = coord_nodos[1]
    tri_x = [p1[0] - sz, p1[0] + sz, p1[0], p1[0] - sz]
    tri_y = [p1[1] - sz*0.8, p1[1] - sz*0.8, p1[1], p1[1] - sz*0.8]
    ax.fill(tri_x, tri_y, color='#e74c3c', alpha=0.7, zorder=3)
    ax.plot([p1[0] - sz*1.3, p1[0] + sz*1.3], [p1[1] - sz*0.8, p1[1] - sz*0.8],
            color='#c0392b', linewidth=2, zorder=3)
    for i in range(5):
        xi = p1[0] - sz*1.3 + i * sz*0.65
        ax.plot([xi, xi - sz*0.3], [p1[1] - sz*0.8, p1[1] - sz*1.3],
                color='#c0392b', linewidth=1, zorder=3)

    p4 = coord_nodos[4]
    tri_x = [p4[0] - sz, p4[0] + sz, p4[0], p4[0] - sz]
    tri_y = [p4[1] - sz*0.8, p4[1] - sz*0.8, p4[1], p4[1] - sz*0.8]
    ax.fill(tri_x, tri_y, color='#2ecc71', alpha=0.7, zorder=3)
    ax.plot([p4[0] - sz*1.3, p4[0] + sz*1.3], [p4[1] - sz*0.8, p4[1] - sz*0.8],
            color='#27ae60', linewidth=2, zorder=3)
    for i in range(5):
        xi = p4[0] - sz*1.3 + i * sz*0.65
        ax.plot([xi, xi - sz*0.3], [p4[1] - sz*0.8, p4[1] - sz*1.3],
                color='#27ae60', linewidth=1, zorder=3)

def dibujar_marco_base(ax, titulo):
    for tag in elementos:
        pI, pJ = obtener_puntos_elemento(tag)
        dibujar_elemento_base(ax, pI, pJ)
    dibujar_nodos(ax)
    dibujar_apoyos(ax)
    ax.set_aspect('equal')
    ax.set_title(titulo, fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel('x (m)', fontsize=9)
    ax.set_ylabel('y (m)', fontsize=9)
    ax.grid(True, alpha=0.2)

def calcular_escala(force_vals, max_px=1.5):
    max_abs = max(abs(v) for v in force_vals)
    if max_abs < 1e-10:
        return 1.0
    return max_px / max_abs

# ============================================================================
# 4. DIAGRAMA DE MOMENTOS FLECTORES (DMF)
# ============================================================================

fig1, ax1 = plt.subplots(figsize=(10, 7))
fig1.patch.set_facecolor('white')

momentos = []
for tag in elementos:
    f = fuerzas[tag]
    momentos.extend([f['Mz_I'], f['Mz_J']])

escala_m = calcular_escala(momentos, max_px=1.8)

dibujar_marco_base(ax1, 'DIAGRAMA DE MOMENTOS FLECTORES (DMF)')

for tag in elementos:
    pI, pJ = obtener_puntos_elemento(tag)
    f = fuerzas[tag]
    dibujar_diag_funcional(ax1, pI, pJ, f['Mz_I'], f['Mz_J'],
                           scale=escala_m, fill_color='#e74c3c', alpha=0.35,
                           line_color='#c0392b')

fig1.text(0.5, 0.01, 'Unidades: tonf*m  |  Convencion: traccion interna positiva',
          ha='center', fontsize=8, style='italic', color='gray')
plt.tight_layout()
plt.subplots_adjust(bottom=0.06)
fig1.savefig('DMF.png', dpi=200, bbox_inches='tight', facecolor='white')
print("Guardado: DMF.png")

# ============================================================================
# 5. DIAGRAMA DE CORTANTES (DEC)
# ============================================================================

fig2, ax2 = plt.subplots(figsize=(10, 7))
fig2.patch.set_facecolor('white')

cortantes = []
for tag in elementos:
    f = fuerzas[tag]
    cortantes.extend([f['Fy_I'], f['Fy_J']])

escala_v = calcular_escala(cortantes, max_px=1.8)

dibujar_marco_base(ax2, 'DIAGRAMA DE ESFUERZOS CORTANTES (DEC)')

for tag in elementos:
    pI, pJ = obtener_puntos_elemento(tag)
    f = fuerzas[tag]
    dibujar_diag_funcional(ax2, pI, pJ, f['Fy_I'], f['Fy_J'],
                           scale=escala_v, fill_color='#3498db', alpha=0.35,
                           line_color='#2980b9')

fig2.text(0.5, 0.01, 'Unidades: tonf  |  Convencion: cortante positivo segun sentido local',
          ha='center', fontsize=8, style='italic', color='gray')
plt.tight_layout()
plt.subplots_adjust(bottom=0.06)
fig2.savefig('DEC.png', dpi=200, bbox_inches='tight', facecolor='white')
print("Guardado: DEC.png")

# ============================================================================
# 6. DIAGRAMA DE AXIALES (DFA)
# ============================================================================

fig3, ax3 = plt.subplots(figsize=(10, 7))
fig3.patch.set_facecolor('white')

axiales = []
for tag in elementos:
    f = fuerzas[tag]
    axiales.extend([f['Fx_I'], f['Fx_J']])

escala_n = calcular_escala(axiales, max_px=1.8)

dibujar_marco_base(ax3, 'DIAGRAMA DE FUERZAS AXIALES (DFA)')

for tag in elementos:
    pI, pJ = obtener_puntos_elemento(tag)
    f = fuerzas[tag]
    dibujar_diag_funcional(ax3, pI, pJ, f['Fx_I'], f['Fx_J'],
                           scale=escala_n, fill_color='#2ecc71', alpha=0.35,
                           line_color='#27ae60')

fig3.text(0.5, 0.01, 'Unidades: tonf  |  Convencion: tension positiva',
          ha='center', fontsize=8, style='italic', color='gray')
plt.tight_layout()
plt.subplots_adjust(bottom=0.06)
fig3.savefig('DFA.png', dpi=200, bbox_inches='tight', facecolor='white')
print("Guardado: DFA.png")

# ============================================================================
# 7. RESUMEN EN CONSOLA
# ============================================================================

print("\n" + "=" * 60)
print("  RESUMEN DE FUERZAS INTERNAS")
print("=" * 60)
for tag, f in fuerzas.items():
    ni, nj = elementos[tag]['nodos']
    print(f"\n  Elemento {tag} ({elementos[tag]['nombre']}):")
    print(f"    Nodo {ni}: Fx={f['Fx_I']:+.4f} | Fy={f['Fy_I']:+.4f} | Mz={f['Mz_I']:+.4f}")
    print(f"    Nodo {nj}: Fx={f['Fx_J']:+.4f} | Fy={f['Fy_J']:+.4f} | Mz={f['Mz_J']:+.4f}")

print("\nDiagramas generados exitosamente!")
