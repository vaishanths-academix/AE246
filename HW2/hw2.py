import math

# Q1(a) Geometry (ft, ft^2, deg)
x_LE_root, x_TE_root, y_root = 134.5380, 148.1480, 6.2679
x_LE_tip,  x_TE_tip,  y_tip  = 147.1589, 152.8730, 24.5466

S_w    = 1950.0     # ft^2 (HW#1)
cbar_w = 16.6417    # ft   (HW#1)
x_cg   = 73.2359    # ft
x_w_ac = 73.5687    # ft

semi_span = y_tip - y_root          # ft
b_h = 2.0 * semi_span               # ft

c_r = x_TE_root - x_LE_root         # ft
c_t = x_TE_tip  - x_LE_tip          # ft
lam = c_t / c_r                     # -

S_h = (c_r + c_t) * semi_span        # ft^2  (since b_h/2 = semi_span)
cbar_h = (2/3) * c_r * (1 + lam + lam**2) / (1 + lam)  # ft

y_mac_local = (semi_span / 3) * (1 + 2*lam) / (1 + lam)  # ft
y_mac_abs = y_root + y_mac_local                          # ft

tan_LE = (x_LE_tip - x_LE_root) / semi_span
Lambda_LE_deg = math.degrees(math.atan(tan_LE))

x_qc_root = x_LE_root + 0.25 * c_r
x_qc_tip  = x_LE_tip  + 0.25 * c_t
tan_qc = (x_qc_tip - x_qc_root) / semi_span
Lambda_qc_deg = math.degrees(math.atan(tan_qc))

x_LE_at_MAC = x_LE_root + (y_mac_abs - y_root) * tan_LE
x_h_ac = x_LE_at_MAC + 0.25 * cbar_h          # ft
l_h = x_h_ac - x_cg                           # ft

V_h = (S_h * l_h) / (S_w * cbar_w)            # -

# Q1(b) Tail lift-curve slope (1/rad)
cl_alpha_2D = 7.2856
CL_alpha_w = 6.4070                           # 1/rad (HW#1)
CL_alpha_h = cl_alpha_2D * math.cos(math.radians(Lambda_qc_deg))  # 1/rad

# Q1(c) Stick-fixed neutral point (ft, -)
deda = 0.5151    # -
eta  = 0.9       # -
Cm_f_alpha = 0.0 # -

delta = (CL_alpha_h / CL_alpha_w) * (1 - deda) * eta * V_h
x_np = x_w_ac + delta * cbar_w
static_margin = (x_np - x_cg) / cbar_w

# Q2
K = 0.0474
K_h = 0.0816
CL = 0.3221

CL_h_del_e = K_h * CL_alpha_h

CL_del_e = eta * (S_h/S_w) * CL_h_del_e

CM_del_e = - eta * V_h * CL_h_del_e

CD_del_e = 2 * K * CL * CL_del_e


# Output
print("==== Q1 Horizontal Tail + Neutral Point (so far) ====\n")

print("Assumptions used:")
print("  • Tail AC at quarter-chord of MAC")
print("  • C_Lα,h ≈ c_lα cos(Λ_c/4)  (same simplification used in HW#1)")
print(f"  • dε/dα = {deda:.4f}")
print(f"  • η = {eta:.3f}")
print(f"  • C_mf,α = {Cm_f_alpha:.1f} (neglect fuselage)\n")

print("Q1(a) Geometry (ft, ft^2, deg):")
print(f"  c_r = x_TE_root - x_LE_root = {x_TE_root:.4f} - {x_LE_root:.4f} = {c_r:.4f} ft")
print(f"  c_t = x_TE_tip  - x_LE_tip  = {x_TE_tip:.4f} - {x_LE_tip:.4f} = {c_t:.4f} ft")
print(f"  λ = c_t/c_r = {lam:.4f}")
print(f"  b_h = 2(y_tip - y_root) = {b_h:.4f} ft")
print(f"  S_h = (c_r + c_t)(b_h/2) = {S_h:.4f} ft^2")
print(f"  c̄_h = {cbar_h:.4f} ft")
print(f"  Λ_LE = {Lambda_LE_deg:.4f} deg")
print(f"  Λ_c/4 = {Lambda_qc_deg:.4f} deg")
print(f"  x_h,ac = {x_h_ac:.4f} ft")
print(f"  l_h = x_h,ac - x_cg = {x_h_ac:.4f} - {x_cg:.4f} = {l_h:.4f} ft")
print(f"  V_h = {V_h:.4f}\n")

print("Q1(b) Lift-curve slopes (1/rad):")
print(f"  C_Lα,w (HW#1) = {CL_alpha_w:.4f}")
print(f"  C_Lα,h = c_lα cos(Λ_c/4) = {cl_alpha_2D:.4f}·cos({Lambda_qc_deg:.4f}°) = {CL_alpha_h:.4f}\n")

print("Q1(c) Stick-fixed neutral point (neglect fuselage):")
print(f"  (1 - dε/dα) = {1 - deda:.4f}")
print(f"  Δ = (C_Lα,h/C_Lα,w)(1-dε/dα)ηV_h = {delta:.4f}")
print(f"  x_np = x_w,ac + Δ·c̄ = {x_w_ac:.4f} + {delta:.4f}·{cbar_w:.4f} = {x_np:.4f} ft")
print(f"  Static margin = (x_np - x_cg)/c̄ = ({x_np:.4f} - {x_cg:.4f})/{cbar_w:.4f} = {static_margin:.4f}")

print("\nStability check:")
print("Statically stable (stick-fixed)" if static_margin > 0 else "  ✘ Not statically stable (stick-fixed)")

print("\n==== Q2 Elevator Control Derivatives (per rad) ====\n")

print(f"CL_h_del_e = {CL_h_del_e:.6f}  # Tail lift derivative wrt δe (1/rad)")
print(f"CL_del_e   = {CL_del_e:.6f}  # Aircraft lift derivative wrt δe (1/rad)")
print(f"CM_del_e   = {CM_del_e:.6f}  # Aircraft pitching moment derivative wrt δe (1/rad)")
print(f"CD_del_e   = {CD_del_e:.6f}  # Aircraft drag derivative wrt δe (1/rad)")