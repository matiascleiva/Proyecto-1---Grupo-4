"""
Diagramas de Esfuerzos - Marco 2D con Asentamientos
====================================================
Estructura: 1 columna (AB) + 1 viga continua (BD)
Punto C es un punto intermedio de la viga, NO una rotula.

Genera:
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

# ============================================================================
# 1. MODELO OPENSEES
# ============================================================================

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 4.0)
ops.node(3, 4.0, 4.0)
ops.node(4, 8.0, 4.0)

ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 0)

A = 1.0e12
E = 1.0
I = 1000.0

ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, A, E, I, 1)
ops.element('elasticBeamColumn', 2, 2, 3, A, E, I, 1)
ops.element('elasticBeamColumn', 3, 3, 4, A, E, I, 1)

ops.constraints('Transformation')
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.sp(1, 2, -0.02)
ops.sp(4, 2, -0.04)

ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-10, 10)
ops.algorithm('Linear')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')
ops.analyze(1)
ops.reactions()

# ============================================================================
# 2. EXTRACCION DE FUERZAS Y AGRUPACION POR MIEMBRO ESTRUCTURAL
# ============================================================================

coord_nodos = {
    1: np.array([0.0, 0.0]),
    2: np.array([0.0, 4.0]),
    3: np.array([4.0, 4.0]),
    4: np.array([8.0, 4.0]),
}

f_elem = {}
for tag in [1, 2, 3]:
    f = ops.eleForce(tag)
    f_elem[tag] = {
        'Fx_I': f[0], 'Fy_I': f[1], 'Mz_I': -f[2],
        'Fx_J': f[3], 'Fy_J': f[4], 'Mz_J': -f[5],
    }

col_AB = {
    'pI': coord_nodos[1],
    'pJ': coord_nodos[2],
    'Fx_I': f_elem[1]['Fx_I'], 'Fy_I': f_elem[1]['Fy_I'], 'Mz_I': f_elem[1]['Mz_I'],
    'Fx_J': f_elem[1]['Fx_J'], 'Fy_J': f_elem[1]['Fy_J'], 'Mz_J': f_elem[1]['Mz_J'],
}

viga_BD = {
    'pI': coord_nodos[2],
    'pC': coord_nodos[3],
    'pJ': coord_nodos[4],
    'Fx_B': f_elem[2]['Fx_I'], 'Fy_B': f_elem[2]['Fy_I'], 'Mz_B': f_elem[2]['Mz_I'],
    'Fx_C': f_elem[2]['Fx_J'],
    'Fy_C': f_elem[2]['Fy_J'],
    'Mz_C': f_elem[2]['Mz_J'],
    'Fx_D': f_elem[3]['Fx_J'], 'Fy_D': f_elem[3]['Fy_J'], 'Mz_D': f_elem[3]['Mz_J'],
}

print("Fuerzas por miembro estructural:")
print(f"  Columna AB: M_A={col_AB['Mz_I']:+.4f} | M_B={col_AB['Mz_J']:+.4f}")
print(f"  Viga BD:    M_B={viga_BD['Mz_B']:+.4f} | M_C={viga_BD['Mz_C']:+.4f} | M_D={viga_BD['Mz_D']:+.4f}")

# ============================================================================
# 3. FUNCIONES DE GRAFICACION
# ============================================================================

def dibujar_elemento(ax, pI, pJ, color='#2c3e50', lw=2.0):
    ax.plot([pI[0], pJ[0]], [pI[1], pJ[1]], color=color, linewidth=lw, zorder=2)

def dibujar_diag(ax, pI, pJ, val_I, val_J, n_pts=50, scale=1.0,
                  fill_color='#3498db', alpha=0.3, line_color='#2980b9',
                  label_val=True):
    dx = pJ[0] - pI[0]
    dy = pJ[1] - pI[1]
    L = np.sqrt(dx**2 + dy**2)
    tx = np.array([dx, dy]) / L
    nx = np.array([-tx[1], tx[0]])

    s = np.linspace(0, 1, n_pts)
    vals = val_I + (val_J - val_I) * s
    pts_elem = np.array([pI + t * (pJ - pI) for t in s])
    pts_diag = pts_elem + np.outer(vals * scale, nx)

    poly_x = [pI[0]] + list(pts_diag[:, 0]) + [pJ[0]]
    poly_y = [pI[1]] + list(pts_diag[:, 1]) + [pJ[1]]
    ax.fill(poly_x, poly_y, color=fill_color, alpha=alpha, zorder=1)
    ax.plot(pts_diag[:, 0], pts_diag[:, 1], color=line_color, linewidth=2, zorder=3)

    if label_val:
        for idx in [0, -1]:
            v = vals[idx]
            pt = pts_diag[idx]
            offset = nx * 0.35 * np.sign(v) if abs(v) > 1e-6 else np.array([0, 0])
            if abs(v) > 1e-6:
                ax.annotate(f'{v:.3f}', xy=(pt[0], pt[1]),
                            xytext=(pt[0] + offset[0], pt[1] + offset[1]),
                            fontsize=8, fontweight='bold', ha='center', va='center',
                            color=line_color,
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                      edgecolor=line_color, alpha=0.8), zorder=5)

def dibujar_diag_3pt(ax, pI, pC, pJ, val_I, val_C, val_J, scale=1.0,
                      fill_color='#3498db', alpha=0.3, line_color='#2980b9',
                      label_val=True):
    dibujar_diag(ax, pI, pC, val_I, val_C, scale=scale,
                 fill_color=fill_color, alpha=alpha, line_color=line_color,
                 label_val=False)
    dibujar_diag(ax, pC, pJ, val_C, val_J, scale=scale,
                 fill_color=fill_color, alpha=alpha, line_color=line_color,
                 label_val=False)
    for pt, val in [(pI, val_I), (pC, val_C), (pJ, val_J)]:
        dx = pJ[0] - pI[0]
        dy = pJ[1] - pI[1]
        L = np.sqrt(dx**2 + dy**2)
        tx = np.array([dx, dy]) / L
        nx = np.array([-tx[1], tx[0]])
        offset = nx * 0.35 * np.sign(val) if abs(val) > 1e-6 else np.array([0, 0])
        if abs(val) > 1e-6:
            pt_diag = pt + val * scale * nx
            ax.annotate(f'{val:.3f}', xy=(pt_diag[0], pt_diag[1]),
                        xytext=(pt_diag[0] + offset[0], pt_diag[1] + offset[1]),
                        fontsize=8, fontweight='bold', ha='center', va='center',
                        color=line_color,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                  edgecolor=line_color, alpha=0.8), zorder=5)

def dibujar_nodos(ax):
    for tag, pos in coord_nodos.items():
        ax.plot(pos[0], pos[1], 'o', color='#2c3e50', markersize=6, zorder=4)
        if tag == 1:
            offset = np.array([-0.25, -0.35])
        elif tag == 4:
            offset = np.array([0.25, -0.35])
        elif tag == 3:
            offset = np.array([0.0, 0.35])
        else:
            offset = np.array([-0.35, 0.0])
        ax.annotate(f'N{tag}', xy=(pos[0], pos[1]),
                     xytext=(pos[0] + offset[0], pos[1] + offset[1]),
                     fontsize=9, fontweight='bold', ha='center', va='center', color='#2c3e50')

def dibujar_punto_C(ax):
    pos = coord_nodos[3]
    ax.plot(pos[0], pos[1], 's', color='#f39c12', markersize=5, zorder=4)
    ax.annotate('C', xy=(pos[0], pos[1]),
                xytext=(pos[0] + 0.3, pos[1] + 0.35),
                fontsize=8, fontstyle='italic', color='#f39c12',
                fontweight='bold', zorder=5)

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
    dibujar_elemento(ax, col_AB['pI'], col_AB['pJ'])
    dibujar_elemento(ax, viga_BD['pI'], viga_BD['pJ'])
    dibujar_nodos(ax)
    dibujar_punto_C(ax)
    dibujar_apoyos(ax)
    ax.set_aspect('equal')
    ax.set_title(titulo, fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel('x (m)', fontsize=9)
    ax.set_ylabel('y (m)', fontsize=9)
    ax.grid(True, alpha=0.2)

def calcular_escala(vals, max_px=1.5):
    max_abs = max(abs(v) for v in vals)
    if max_abs < 1e-10:
        return 1.0
    return max_px / max_abs

# ============================================================================
# 4. DIAGRAMA DE MOMENTOS FLECTORES (DMF)
# ============================================================================

fig1, ax1 = plt.subplots(figsize=(10, 7))
fig1.patch.set_facecolor('white')

m_vals = [col_AB['Mz_I'], col_AB['Mz_J'],
          viga_BD['Mz_B'], viga_BD['Mz_C'], viga_BD['Mz_D']]
escala_m = calcular_escala(m_vals, max_px=1.8)

dibujar_marco_base(ax1, 'DIAGRAMA DE MOMENTOS FLECTORES (DMF)')

dibujar_diag(ax1, col_AB['pI'], col_AB['pJ'],
             col_AB['Mz_I'], col_AB['Mz_J'], scale=escala_m,
             fill_color='#e74c3c', alpha=0.35, line_color='#c0392b')

dibujar_diag_3pt(ax1, viga_BD['pI'], viga_BD['pC'], viga_BD['pJ'],
                 viga_BD['Mz_B'], viga_BD['Mz_C'], viga_BD['Mz_D'],
                 scale=escala_m, fill_color='#e74c3c', alpha=0.35, line_color='#c0392b')

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

v_vals = [col_AB['Fy_I'], col_AB['Fy_J'],
          viga_BD['Fy_B'], viga_BD['Fy_C'], viga_BD['Fy_D']]
escala_v = calcular_escala(v_vals, max_px=1.8)

dibujar_marco_base(ax2, 'DIAGRAMA DE ESFUERZOS CORTANTES (DEC)')

dibujar_diag(ax2, col_AB['pI'], col_AB['pJ'],
             col_AB['Fy_I'], col_AB['Fy_J'], scale=escala_v,
             fill_color='#3498db', alpha=0.35, line_color='#2980b9')

dibujar_diag_3pt(ax2, viga_BD['pI'], viga_BD['pC'], viga_BD['pJ'],
                 viga_BD['Fy_B'], viga_BD['Fy_C'], viga_BD['Fy_D'],
                 scale=escala_v, fill_color='#3498db', alpha=0.35, line_color='#2980b9')

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

n_vals = [col_AB['Fx_I'], col_AB['Fx_J'],
          viga_BD['Fx_B'], viga_BD['Fx_C'], viga_BD['Fx_D']]
escala_n = calcular_escala(n_vals, max_px=1.8)

dibujar_marco_base(ax3, 'DIAGRAMA DE FUERZAS AXIALES (DFA)')

dibujar_diag(ax3, col_AB['pI'], col_AB['pJ'],
             col_AB['Fx_I'], col_AB['Fx_J'], scale=escala_n,
             fill_color='#2ecc71', alpha=0.35, line_color='#27ae60')

dibujar_diag_3pt(ax3, viga_BD['pI'], viga_BD['pC'], viga_BD['pJ'],
                 viga_BD['Fx_B'], viga_BD['Fx_C'], viga_BD['Fx_D'],
                 scale=escala_n, fill_color='#2ecc71', alpha=0.35, line_color='#27ae60')

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
print("  RESUMEN DE FUERZAS POR MIEMBRO ESTRUCTURAL")
print("=" * 60)

print(f"\n  Columna AB:")
print(f"    A: Fx={col_AB['Fx_I']:+.4f} | Fy={col_AB['Fy_I']:+.4f} | Mz={col_AB['Mz_I']:+.4f}")
print(f"    B: Fx={col_AB['Fx_J']:+.4f} | Fy={col_AB['Fy_J']:+.4f} | Mz={col_AB['Mz_J']:+.4f}")

print(f"\n  Viga BD (continua, sin rotula en C):")
print(f"    B: Fx={viga_BD['Fx_B']:+.4f} | Fy={viga_BD['Fy_B']:+.4f} | Mz={viga_BD['Mz_B']:+.4f}")
print(f"    C: Fx={viga_BD['Fx_C']:+.4f} | Fy={viga_BD['Fy_C']:+.4f} | Mz={viga_BD['Mz_C']:+.4f}")
print(f"    D: Fx={viga_BD['Fx_D']:+.4f} | Fy={viga_BD['Fy_D']:+.4f} | Mz={viga_BD['Mz_D']:+.4f}")

print("\nDiagramas generados exitosamente!")
