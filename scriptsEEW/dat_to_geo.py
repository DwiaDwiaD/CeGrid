import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev

# WORKS ONLY FOR EPPLER61 RIGHT NOW!!!!!

# -------------------------------------------------
# Controls
# -------------------------------------------------

if len(sys.argv) < 5:
    print("Usage (Uniform):  python3 dat_to_geo.py {filename} {pathToWD} {BL_UNIFORM}")
    print("Usage (Variable): python3 dat_to_geo.py {filename} {pathToWD} {BL_TE_UPPER} {BL_TE_LOWER}")
    sys.exit(1)

name = sys.argv[1]
SCRIPT_DIR = sys.argv[2]

# Flexible arguments to support both uniform and variable thickness
BL_TE_UPPER = float(sys.argv[3])
GROUND = float(sys.argv[4])
# GNDClr = 0.01
GNDClr = 0
AOA_DEG = float(sys.argv[5])
first_layerH = float(sys.argv[6])
NUMlayers = float(sys.argv[7])
PLOT = float(sys.argv[8])
NRESAMPLE = int(sys.argv[9])

GROWTH_EXP = 0.1   # 1.0 = linear BL thickness growth, < 1.0 = rapid initial growth
te_point = (1.0, 0.0) # assumed trailing edge is the default (1.0,0.0) in the airfoil coordinates!

BL_TE_LOWER = min(GROUND-GNDClr,BL_TE_UPPER) 
BL_LE = min(BL_TE_LOWER + 0.9 * np.sin(np.radians(AOA_DEG)), BL_TE_UPPER)

# -------------------------------------------------
# Functions
# -------------------------------------------------
from scipy.optimize import brentq

def find_growth_rate(H, y1, N):
    # Objective function: Total Thickness Error = 0
    def objective(r):
        if abs(r - 1.0) < 1e-9:
            return y1 * N - H
        return y1 * (1 - r**N) / (1 - r) - H
    
    # Growth rate is usually between 1.01 and 1.5
    return brentq(objective, 0.95, 2.0)


def bunching_array(n, power_start, power_end):
    """
    Maps a linear space to a clustered space.
    High power_start clusters near the start index. 
    High power_end clusters near the end index.
    """
    t = np.linspace(0, 1, n)
    res = np.zeros(n)
    for i in range(1, n - 1):
        res[i] = (t[i]**power_start) / (t[i]**power_start + (1.0 - t[i])**power_end)
    res[-1] = 1.0
    return res

def Translate(translate, coords):
    tx, ty = translate
    xs, ys = coords
    xs = np.array(xs)
    ys = np.array(ys)
    xs += tx
    ys += ty
    coords = (xs, ys)
    return coords

def Rotate(angle, coords, point=(0, 0)):
    # Move rotation point to origin
    coords = Translate((-point[0], -point[1]), coords)

    xs, ys = coords
    angle_rad = np.radians(angle)

    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)

    # Rotate about origin
    xs_rot = xs * cos_a - ys * sin_a
    ys_rot = xs * sin_a + ys * cos_a

    coords = (xs_rot, ys_rot)

    # Move back
    coords = Translate(point, coords)

    return coords

def evalNormals(spline, coords):
    tck, u_spline = spline
    x_spline, y_spline = coords
    dx, dy = splev(u_spline, tck, der=1)
    mag = np.hypot(dx, dy)
    nx = -dy / mag
    ny = dx / mag

    # Ensure normals point outward
    cx, cy = np.mean(x_spline), np.mean(y_spline)
    dot_products = nx * (x_spline - cx) + ny * (y_spline - cy)
    if np.sum(dot_products) < 0:
        nx, ny = -nx, -ny
    
    return (nx,ny)

def ArcLens(coords): #returns array of cumulative lengths
    xs, ys = coords
    n = len(xs)
    arclength = np.zeros(n)
    
    for i in range(1,n):
        length = np.hypot(xs[i] - xs[i-1], ys[i] - ys[i-1])
        arclength[i] = arclength[i-1] + length
    arclength /= np.max(arclength)
    return arclength

# -------------------------------------------------
# Growth Rate
# -------------------------------------------------

# Example usage:
# H = BL_TE_UPPER (total thickness), y1 = first_layerH (first layer), N = 100 (layers)
growthR_TEup = find_growth_rate(H=BL_TE_UPPER, y1=first_layerH, N=NUMlayers)
growthR_TElow = find_growth_rate(H=BL_TE_LOWER, y1=first_layerH, N=NUMlayers)
growthR_LE = find_growth_rate(H=BL_LE, y1=first_layerH, N=NUMlayers)


# -------------------------------------------------
# Read Data & Spline
# -------------------------------------------------
raw = np.loadtxt(name)
if np.allclose(raw[0], raw[-1]):
    raw = raw[:-1]


# splprep fits a smooth curve parameterized by u (0 to 1)
tck, u_raw = splprep([raw[:, 0], raw[:, 1]], s=0.0, k=3, per=False)

u_new = u_raw
x_new, y_new = splev(u_new, tck)

# -------------------------------------------------
# Strictly Enforce (1.0, 0.0) at Trailing Edges
# -------------------------------------------------
x_new[0], y_new[0] = te_point[0], te_point[1]
x_new[-1], y_new[-1] = te_point[0], te_point[1]

thick_hi = int(np.argmax(y_new))
thick_lo = int(np.argmin(y_new))

u_thick_hi = float(u_new[thick_hi])
u_thick_lo = float(u_new[thick_lo])


if thick_hi<thick_lo: #anticlockwise, TE-LE-TE
    x_top = x_new[:thick_hi]
    y_top = y_new[:thick_hi]
    
    x_bot = x_new[thick_lo:]
    y_bot = y_new[thick_lo:]  

if thick_hi>thick_lo: #clockwise, TE-LE-TE
    x_top = x_new[thick_hi:]
    y_top = y_new[thick_hi:]
    
    x_bot = x_new[:thick_lo]
    y_bot = y_new[:thick_lo]

# tck_top, u_top = splprep([x_top, y_top], s=0.0, k=3, per=False)
# tck_bot, u_bot = splprep([x_bot, y_bot], s=0.0, k=3, per=False)

bot_angle = AOA_DEG - 0.5 # this angle should be such that the leading end does not strike the ground(perhaps use BL_LE to figure out??)
x_off_bot, y_off_bot = Translate((0,-BL_TE_LOWER),Rotate(bot_angle, (x_bot, y_bot), (1,0)))

x_off_top, y_off_top = Translate((0,BL_TE_UPPER),Rotate(bot_angle, (x_top, y_top), (1,0)))

if (x_top[0] == 1) and (y_top[0] == 0): # i.e. if anticlockwise, TE-LE-TE
    BL_hi = np.hypot((x_top[-1]-x_off_top[-1]),(y_top[-1]-y_off_top[-1]))
    BL_lo = np.hypot((x_bot[0]-x_off_bot[0]),(y_top[0]-y_off_bot[0]))

else: # i.e. if clockwise, TE-LE-TE
    BL_hi = np.hypot((x_top[0]-x_off_top[0]),(y_top[0]-y_off_top[0]))
    BL_lo = np.hypot((x_bot[-1]-x_off_bot[-1]),(y_top[-1]-y_off_bot[-1]))


x_nose = x_new[min(thick_lo, thick_hi):max(thick_lo, thick_hi)]
y_nose = y_new[min(thick_lo, thick_hi):max(thick_lo, thick_hi)]

x_nose, y_nose = Rotate(-AOA_DEG, (x_nose, y_nose)) #this rotate so that we can check ground level using just y
tck_nose, u_nose = splprep([x_nose, y_nose], s=0.0, k=3, per=False)
nx_nose, ny_nose = evalNormals((tck_nose, u_nose),(x_nose, y_nose))

if (x_top[0] == 1) and (y_top[0] == 0): # i.e. if anticlockwise, TE-LE-TE
    x_nose,y_nose = (x_nose[::-1],y_nose[::-1])
    nx_nose, ny_nose = (nx_nose[::-1],ny_nose[::-1])

N_nose = len(x_nose)

arclens = ArcLens((x_nose,y_nose))

# print(arclens)

# -------------------------------------------------
# Boundary-layer thickness over the nose
# -------------------------------------------------

Dists = np.zeros(N_nose)

for i in range(N_nose):

    # arclens already runs from 0 -> 1
    s = arclens[i]

    # desired thickness between lower and upper values
    target = BL_lo + (BL_hi - BL_lo) * s**GROWTH_EXP

    # maximum distance before the offset would touch the ground
    if ny_nose[i] < 0:
        ground_limit = (GROUND + y_nose[i]) / (-ny_nose[i])
    else:
        ground_limit = 1e9

    Dists[i] = min(target, ground_limit)

# -------------------------------------------------
# Smooth the distance field
# -------------------------------------------------

for _ in range(5):
    Dists[1:-1] = (
        0.25*Dists[:-2]
        + 0.50*Dists[1:-1]
        + 0.25*Dists[2:]
    )

# keep end values fixed
Dists[0]  = BL_lo
Dists[-1] = BL_hi

# -------------------------------------------------
# Offset nose
# -------------------------------------------------

x_nose = x_nose + Dists*nx_nose
y_nose = y_nose + Dists*ny_nose

if (x_top[0] == 1) and (y_top[0] == 0): # i.e. if anticlockwise, TE-LE-TE
    x_nose = x_nose[::-1]
    y_nose = y_nose[::-1]

x_off_nose, y_off_nose = Rotate(AOA_DEG, (x_nose, y_nose))

translate = (-te_point[0],-te_point[1])
x_off_nose, y_off_nose = Translate(translate,(x_off_nose,y_off_nose))
x_off_top, y_off_top = Translate(translate,(x_off_top,y_off_top))
x_off_bot, y_off_bot = Translate(translate,(x_off_bot,y_off_bot))

x_off_nose, y_off_nose = Rotate(-AOA_DEG, (x_off_nose, y_off_nose))
x_off_top,y_off_top= Rotate(-AOA_DEG, (x_off_top,y_off_top))
x_off_bot,y_off_bot= Rotate(-AOA_DEG, (x_off_bot,y_off_bot))

y_off_bot = -BL_TE_LOWER*np.ones_like(y_off_bot)

if (x_top[0] == 1) and (y_top[0] == 0): # i.e. if anticlockwise, TE-LE-TE
    x_off = np.concatenate((x_off_top, x_off_nose, x_off_bot))
    y_off = np.concatenate((y_off_top, y_off_nose, y_off_bot))

else: # i.e. if clockwise, TE-LE-TE
    x_off = np.concatenate((x_off_bot, x_off_nose, x_off_top))
    y_off = np.concatenate((y_off_bot, y_off_nose, y_off_top))

tck_off, _ = splprep([x_off, y_off], s=0.0, k=3, per=False)
u_off = np.linspace(0.0, 1.0, NRESAMPLE)
x_off, y_off = splev(u_off, tck_off)


# -------------------------------------------------
# Resample with Aggressive LE Bunching
# -------------------------------------------------
u_dense = np.linspace(0.0, 1.0, 10000)
x_dense, _ = splev(u_dense, tck)
u_le = float(u_dense[np.argmin(x_dense)])

LE_BUNCHING = 2.0  # Higher = significantly denser at the Leading Edge
TE_BUNCHING = 1.0  # Higher = denser at the Trailing Edge

n_up = NRESAMPLE // 2
n_low = NRESAMPLE - n_up

# Upper Surface: TE (u=0) to LE (u=u_le)
u_up_frac = bunching_array(n_up, TE_BUNCHING, LE_BUNCHING)
u_up = u_up_frac * u_le

# Lower Surface: LE (u=u_le) to TE (u=1)
u_low_frac = bunching_array(n_low, LE_BUNCHING, TE_BUNCHING)
u_low = u_le + u_low_frac * (1.0 - u_le)

# Combine (avoiding duplicating the LE point)
u_new = np.concatenate([u_up[:-1], u_low])
x_new, y_new = splev(u_new, tck)

n = len(x_new)

# -------------------------------------------------
# Split by arc length instead of LE index
# -------------------------------------------------

arc_air = ArcLens((x_new, y_new))
arc_off = ArcLens((x_off, y_off))

# Note: Gmsh uses a negative angle in the script, which rotates clockwise.
# Rotate Airfoil Points to Global Frame
translate = (-te_point[0],-te_point[1])
x_new, y_new = Translate(translate,(x_new,y_new))
# REENFORCING TE
x_new[0],y_new[0] = (0,0)

# translate = (-te_point[0],-te_point[1])
# x_off, y_off = Translate(translate,(x_off,y_off))
x_glob, y_glob = Rotate(-AOA_DEG, (x_new,y_new))

# Rotate Boundary Layer Edge Points to Global Frame
# x_off_glob, y_off_glob = Rotate(-AOA_DEG, (x_off,y_off))
x_off_glob, y_off_glob = (x_off,y_off)
# y_off_glob[-1, -len(y_bot)]
split_air = n_up
split_off = np.argmin(np.abs(y_off_glob - np.sin(np.radians(AOA_DEG))))
N_UP = int(NRESAMPLE*arc_off[split_off])
N_UP = NRESAMPLE//2
N_LOW = NRESAMPLE - N_UP

ground_pts = np.where(
    y_off <= -BL_TE_LOWER * (1 - 10**-3)
)[0]

i_left = ground_pts[np.argmin(x_off[ground_pts])]

# Clamp the remaining lower offset boundary to the ground
y_off[i_left:] = -BL_TE_LOWER
y_off_glob[i_left:] = -BL_TE_LOWER

gnd_pt = n + 1 + i_left

# -------------------------------------------------
# Distribute lower BL nodes between curves 4 and 5
# -------------------------------------------------

s_LE  = arc_off[split_off]
s_GND = arc_off[i_left]

# Normalized arc lengths of curves 4 and 5
L4 = s_GND - s_LE
L5 = 1.0 - s_GND

Ltotal = L4 + L5

# N_LOW is the desired TOTAL number of nodes
# along curves 4 + 5, including the shared GND point.
N_segments = N_LOW - 1

# Allocate elements proportionally to physical arc length
Nseg4 = int(round(N_segments * L4 / Ltotal))
Nseg5 = N_segments - Nseg4

# Convert element counts to node counts
N4 = Nseg4 + 1
N5 = Nseg5 + 1

# Airfoil
airfoil_up_ids = list(range(1, split_air + 2))
airfoil_low_ids = list(range(split_air + 1, n)) + [1]

# Offset
off_up_ids = list(range(n + 1, n + split_off + 2))
off_low_ids = list(range(n + split_off + 1,gnd_pt + 1))
off_gnd_ids = list(range(gnd_pt, 2*n + 2))

with open(f"{SCRIPT_DIR}/scriptsEEW/Airfoil_points.geo", "w") as f:
    f.write("// Airfoil points\n")
    for i, (xv, yv) in enumerate(zip(x_glob, y_glob)):
        f.write(f"Point({i + 1}) = {{{xv:.6f}, {yv:.6f}, 0, 1.0}};\n")

    f.write("\n// Offset boundary points\n")
    for i, (xv, yv) in enumerate(zip(x_off, y_off)):
        f.write(f"Point({n + i + 1}) = {{{xv:.6f}, {yv:.6f}, 0, 1.0}};\n")

    f.write(f"\nBSpline(1) = {{{','.join(map(str, airfoil_up_ids))}}}; // Airfoil Upper\n")
    f.write(f"BSpline(2) = {{{','.join(map(str, airfoil_low_ids))}}}; // Airfoil Lower\n")
    f.write(f"BSpline(3) = {{{','.join(map(str, off_up_ids))}}}; // Offset Upper\n")
    f.write(f"BSpline(4) = {{{','.join(map(str, off_low_ids))}}}; // Offset Lower\n")
    f.write(f"BSpline(5) = {{{','.join(map(str, off_gnd_ids))}}}; // Offset Ground\n")

    # Export dynamic point counts so Transfinite Meshing matches exactly
    f.write(f"\nN_UP = {N_UP};\n")
    f.write(f"N_LOW = {N_LOW};\n")
    # f.write(f"\nN_UP = {int(n_up)};\n")
    # f.write(f"N_LOW = {int(n_low)};\n")

    f.write(f"\nLEpoint = {split_air + 1};\n")
    f.write(f"LEoff = {n + split_off + 1};\n")
    f.write(f"TEpoint = {1};\n")
    f.write(f"TEoff_up = {n + 1};\n")
    f.write(f"TEoff_low = {2 * n + 1};\n")
    f.write(f"GNDpoint = {gnd_pt};\n")
    f.write(f"BLAirfoilUp = {BL_TE_UPPER:.4};\n")
    f.write(f"BLAirfoilLow = {BL_TE_LOWER:.4};\n")
    f.write(f"hc = {GROUND};\n")
    f.write(f"grTEup = {growthR_TEup:.6f};\n")
    f.write(f"grTElow = {growthR_TElow:.6f};\n")
    f.write(f"grLE = {growthR_LE:.6f};\n")
    f.write(f"NUMlayers = {NUMlayers};\n")
    f.write(f"N4={N4}; \nN5={N5};\n")

# -------------------------------------------------
# Diagnostic Plotting (Global Frame Rendering)
# -------------------------------------------------
if PLOT:
    # 2D Rotation Matrix to match Gmsh's Rotate {{0,0,1}, {0,0,0}, -AoA_rad}
    aoa_rad = np.radians(AOA_DEG)

    # --- Plotting ---
    plt.figure(figsize=(10, 5))
    plt.plot(x_glob, y_glob, 'k.-', linewidth=1.5, markersize=3, label='Resampled Airfoil (Global)')
    plt.plot(x_off_glob, y_off_glob, 'r.-', linewidth=1.5, markersize=3, label='Boundary Layer Edge (Global)')
    plt.plot(x_off_glob[i_left],y_off_glob[i_left], 's')

    # Draw tie-lines cleanly in the Global Frame
    for i in range(0, n, 4):
        plt.plot([x_glob[i], x_off_glob[i]], [y_glob[i], y_off_glob[i]], color='gray', alpha=0.4, linewidth=0.8)

    # --- DRAW THE GROUND LINE (Dead simple in Global Frame) ---
    x_span = np.linspace(np.min(x_off_glob) - 0.2, np.max(x_off_glob) + 0.5, 100)
    y_ground = np.full_like(x_span, -GROUND)  # Horizontal flat line

    plt.plot(x_span, y_ground, color='brown', linestyle='--', linewidth=2.0, label=f'Ground Plane (y = {-GROUND})')
    # -----------------------------------------------------------

    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.xlabel('Global X')
    plt.ylabel('Global Y')
    plt.title(f'Global Frame Check (AoA = {AOA_DEG}°, Ground H/C = {GROUND})')
    plt.tight_layout()
    plt.show()