"""
OpenSeesPy Script: Marco Hiperestatico 2D con Asentamientos
============================================================
Analisis de un marco 2D sometido a descensos en sus apoyos.
El marco consta de una columna y dos segmentos de viga.

Unidades:
- Fuerzas: tonf (tonelada-fuerza)
- Longitudes: m (metros)
"""

import openseespy.opensees as ops

# ============================================================================
# 1. INICIALIZACION
# ============================================================================

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# ============================================================================
# 2. DEFINICION DE NODOS
# ============================================================================

ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 4.0)
ops.node(3, 4.0, 4.0)
ops.node(4, 8.0, 4.0)

# ============================================================================
# 3. CONDICIONES DE BORDE (APOYOS)
# ============================================================================

ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 0)

# ============================================================================
# 4. PROPIEDADES DE SECCION Y MATERIAL
# ============================================================================

A = 1.0e12
E = 1.0
I = 1000.0

# ============================================================================
# 5. TRANSFORMACION GEOMETRICA
# ============================================================================

ops.geomTransf('Linear', 1)

# ============================================================================
# 6. DEFINICION DE ELEMENTOS
# ============================================================================

ops.element('elasticBeamColumn', 1, 1, 2, A, E, I, 1)
ops.element('elasticBeamColumn', 2, 2, 3, A, E, I, 1)
ops.element('elasticBeamColumn', 3, 3, 4, A, E, I, 1)

# ============================================================================
# 7. ASENTAMIENTOS (DESPLAZAMIENTOS IMPUESTOS)
# ============================================================================

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

# ============================================================================
# 9. EXTRACCION DE RESULTADOS
# ============================================================================

print("\n" + "=" * 60)
print("  RESULTADOS: Marco 2D con Asentamientos")
print("=" * 60)

# --- 9.1 Reacciones ---
print("\n--- REACCIONES ---")
ops.reactions()

node1_rx = ops.nodeReaction(1, 1)
node1_ry = ops.nodeReaction(1, 2)
node1_mz = -ops.nodeReaction(1, 3)

node4_rx = ops.nodeReaction(4, 1)
node4_ry = ops.nodeReaction(4, 2)
node4_mz = -ops.nodeReaction(4, 3)

print(f"\nNodo 1 (Apoyo A - Empotrado):")
print(f"  Fx = {node1_rx:12.6f} tonf")
print(f"  Fy = {node1_ry:12.6f} tonf")
print(f"  Mz = {node1_mz:12.6f} tonf*m")

print(f"\nNodo 4 (Apoyo D - Articulado):")
print(f"  Fx = {node4_rx:12.6f} tonf")
print(f"  Fy = {node4_ry:12.6f} tonf")
print(f"  Mz = {node4_mz:12.6f} tonf*m")

print(f"\nVerificacion de equilibrio:")
print(f"  Sum Fx = {node1_rx + node4_rx:12.6f} tonf (debe ser 0)")
print(f"  Sum Fy = {node1_ry + node4_ry:12.6f} tonf (debe ser 0)")

# --- 9.2 Fuerzas Internas (Momentos Flectores) ---
print("\n--- MOMENTOS FLECTORES ---")
for ele_tag in range(1, 4):
    forces = ops.eleForce(ele_tag)
    print(f"Elemento {ele_tag}: M_I = {-forces[2]:12.6f} tonf*m | M_J = {-forces[5]:12.6f} tonf*m")

# --- 9.3 Desplazamientos ---
print("\n--- DESPLAZAMIENTOS ---")
disp_C = ops.nodeDisp(3, 2)
print(f"Nodo 3 (Punto C): dVc = {disp_C:12.6f} m ({disp_C * 1000:8.3f} mm)")

print("\nDesplazamientos completos por nodo:")
for node_tag in range(1, 5):
    dx = ops.nodeDisp(node_tag, 1)
    dy = ops.nodeDisp(node_tag, 2)
    rz = ops.nodeDisp(node_tag, 3)
    print(f"  Nodo {node_tag}: dx = {dx:12.6f} m | dy = {dy:12.6f} m | rz = {rz:12.6f} rad")

print("\n" + "=" * 60)
print("  Analisis completado exitosamente!")
print("=" * 60)
