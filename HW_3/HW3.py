import math

# =============================================================================
# GIVEN VALUES
# =============================================================================

# Vertical stabiliser geometry points (x, y) in inches
LE_root = (123.3152, 10.6970)
LE_tip  = (147.4552, 32.9545)
TE_root = (146.6092, 10.6970)
TE_tip  = (155.2385, 32.9545)

# Aircraft-level parameters
x_ac_wing     = 73.5687   # ft
x_cg_aircraft = 73.2359   # ft
S_wing        = 1950       # ft^2
b_wing        = 124.2572   # ft
cl_beta       = 7.2856
eta           = 0.9992
CL_alpha_wing = 6.4070
K = 0.0474
CL_wing = 0.32284649


# Horizontal stabiliser parameters
CL_alpha_ht = 6.2955
b_ht        = 36.5574
lambda_ht   = 0.4198
S_h         = 353.2194

# Wing geometry for Q3
wing_LE_root = (81.3008, 46.7028)
wing_LE_tip  = (89.6503, 62.1286)
wing_TE_root = (91.5173, 46.7028)
wing_TE_tip  = (95.3544, 62.1286)

aileron_root_LE = (89.2962, 46.7028)
aileron_tip_LE  = (94.1339, 62.1286)

wing_dihedral = math.radians(5)
ht_dihedral   = math.radians(7)

# Rudder parameters (Q2)
tau = 0.55
K_v = 0.1078

# Sideslip cases (Q2)
beta_1 = 0
beta_2 = math.radians(2.0)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def chord(LE, TE):
    return TE[0] - LE[0]

def span(root, tip):
    return tip[1] - root[1]

def trapezoid_area(c_root, c_tip, b):
    return 0.5 * (c_root + c_tip) * b

def taper_ratio(c_tip, c_root):
    return c_tip / c_root

def mean_aerodynamic_chord(c_root, lam):
    return (2/3) * c_root * (lam**2 + lam + 1) / (lam + 1)

def quarter_chord_x(LE, c):
    return LE[0] + 0.25 * c

def sweep_angle_rad(x_qc_root, x_qc_tip, b):
    return math.atan((x_qc_tip - x_qc_root) / b)

# =============================================================================
# Q1 — VERTICAL STABILISER GEOMETRY AND AERODYNAMICS
# =============================================================================

c_root    = chord(LE_root, TE_root)
c_tip     = chord(LE_tip,  TE_tip)

b_v       = span(LE_root, LE_tip)
S_v       = trapezoid_area(c_root, c_tip, b_v)
lam       = taper_ratio(c_tip, c_root)
c_bar_v   = mean_aerodynamic_chord(c_root, lam)

x_qc_root = quarter_chord_x(LE_root, c_root)
x_qc_tip  = quarter_chord_x(LE_tip,  c_tip)
sweep_rad = sweep_angle_rad(x_qc_root, x_qc_tip, b_v)
sweep_deg = math.degrees(sweep_rad)

x_ac_v    = (((c_bar_v - c_tip) / (c_root - c_tip)) * (x_qc_root - x_qc_tip)) + x_qc_tip
l_v       = x_ac_v - x_cg_aircraft
V_v       = (l_v * S_v) / (S_wing * b_wing)

# Part (b)
CY_v_beta = cl_beta * math.cos(sweep_rad)
C_n_beta  = CY_v_beta * V_v

# =============================================================================
# Q2 — RUDDER DERIVATIVES
# =============================================================================

CL_v_del_r      = K_v * CY_v_beta
do_CY_by_do_delr = CY_v_beta * tau

Sv_Sw    = S_v / S_wing
CY_del_r = Sv_Sw * do_CY_by_do_delr
Cn_del_r = -do_CY_by_do_delr * V_v

def CD_del_r(beta: float) -> float:
    return (2.0 * K_v * CY_v_beta * beta * do_CY_by_do_delr) * Sv_Sw

CD_del_r_case_1 = CD_del_r(beta_1)
CD_del_r_case_2 = CD_del_r(beta_2)
CD_del_r_func = (2.0 * K_v * CY_v_beta * do_CY_by_do_delr) * Sv_Sw

# =============================================================================
# Q3 — WING / HT LATERAL STABILITY (DIHEDRAL EFFECT)
# =============================================================================

#c_root_wing_section = chord(wing_LE_root, wing_TE_root)

c_tip_wing          = chord(wing_LE_tip,  wing_TE_tip)
c_root_wing         = (2 * S_wing / b_wing) - c_tip_wing
lambda_wing         = c_tip_wing / c_root_wing

y_w_ac  = ((b_wing / 2) / 3) * ((1 + 2 * lambda_wing) / (1 + lambda_wing))
y_ht_ac = ((b_ht   / 2) / 3) * ((1 + 2 * lambda_ht)   / (1 + lambda_ht))

Cl_beta_wing = -CL_alpha_wing * wing_dihedral * y_w_ac / b_wing
Cl_beta_ht   = -CL_alpha_ht   * ht_dihedral   * eta * (S_h * y_ht_ac) / (S_wing * b_wing)
Cl_beta      =  Cl_beta_wing + Cl_beta_ht

c_root_aileron = chord(aileron_root_LE, wing_TE_root)
c_tip_aileron  = chord(aileron_tip_LE,  wing_TE_tip)
b_aileron      = span(aileron_root_LE,  aileron_tip_LE)
S_w_aileron    = (c_root_aileron + c_tip_aileron) * b_aileron

y_a = 55.17 #ft
tau_aileron = 0.55

do_CLw_by_do_del_a = CL_alpha_wing * tau_aileron

geo_factor = (y_a * S_w_aileron) / (S_wing * b_wing)

Cl_del_a = do_CLw_by_do_del_a * geo_factor

#Cn_del_a
do_CD_w_by_do_del_a = 2 * K * CL_wing * do_CLw_by_do_del_a
Cn_del_a = - do_CD_w_by_do_del_a * geo_factor

# =============================================================================
# RESULTS
# =============================================================================

print("=" * 45)
print("Q1 - VERTICAL STABILISER")
print("=" * 45)
print(f"  Root chord:        {c_root:.4f} ft")
print(f"  Tip chord:         {c_tip:.4f} ft")
print(f"  Span (b_v):        {b_v:.4f} ft")
print(f"  Area (S_v):        {S_v:.4f} ft^2")
print(f"  Taper ratio:       {lam:.4f}")
print(f"  MAC (c_bar_v):     {c_bar_v:.4f} ft")
print(f"  QC x root:         {x_qc_root:.4f} ft")
print(f"  QC x tip:          {x_qc_tip:.4f} ft")
print(f"  Sweep (c/4):       {sweep_deg:.4f} deg")
print(f"  l_v:               {l_v:.4f} ft")
print(f"  V_v:               {V_v:.4f}")
print(f"  CY_v_beta:         {CY_v_beta:.4f}")
print(f"  Cn_v_beta:         {C_n_beta:.4f}")
print()
print("=" * 45)
print("Q2 - RUDDER DERIVATIVES")
print("=" * 45)
print(f"  CY_del_r:          {CY_del_r:.4f}")
print(f"  Cn_del_r:          {Cn_del_r:.4f}")
print(f"  CD_del_r (beta=0 deg):  {CD_del_r_case_1:.4f}")
print(f"  CD_del_r (beta=2 deg):  {CD_del_r_case_2:.4f}")
print(f"  CD_del_r as a function of beta:  {CD_del_r_func:.4f} * beta")
print()
print("=" * 45)
print("Q3 - LATERAL STABILITY (DIHEDRAL EFFECT)")
print("=" * 45)
#print(f"  y_w_ac:            {y_w_ac:.4f} ft")
#print(f"  y_ht_ac:           {y_ht_ac:.4f} ft")
#print(f"  Cl_beta_wing:      {Cl_beta_wing:.4f}")
#print(f"  Cl_beta_ht:        {Cl_beta_ht:.4f}")
print(f"  Cl_beta:           {Cl_beta:.4f}")
#print(f"  b_aileron:         {b_aileron:.4f} ft")
#print(f"  S_w_aileron:       {S_w_aileron:.4f} ft^2")
print(f"  Cl_delta_a:        {Cl_del_a:.4f}")
print(f"  Cn_delta_a:        {Cn_del_a:.4f}")