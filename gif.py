"""
make_gif.py
Reads frames.bin produced by sheet.c and generates a 3D surface GIF.

Usage:
    pip install matplotlib numpy pillow
    python make_gif.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import struct

# ------------------------------------------------------------------ #
#  Parameters (must match sheet.c)                                     #
# ------------------------------------------------------------------ #
BIN_FILE   = "frames.bin"
OUTPUT_GIF = "output.gif"
L          = 1    # side length [m]
Z_MIN      = -0.1 # fixed z axis min [m]  (adjust if AMP changes)
Z_MAX      =  0.1 # fixed z axis max [m]

# ------------------------------------------------------------------ #
#  Read binary file                                                    #
# ------------------------------------------------------------------ #
with open(BIN_FILE, "rb") as f:
    nframes_header, N = struct.unpack("ii", f.read(8))
    raw = np.frombuffer(f.read(), dtype=np.float32)

# compute actual number of complete frames in the file
nframes = len(raw) // (N * N)
print(f"Header said {nframes_header} frames, got {nframes} complete frames, N={N}")

frames = raw[:nframes * N * N].reshape((nframes, N, N))

# ------------------------------------------------------------------ #
#  Build coordinate grids                                              #
# ------------------------------------------------------------------ #
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)   # both shape (N, N), X[j,i]=x_i, Y[j,i]=y_j

# ------------------------------------------------------------------ #
#  Set up figure                                                       #
# ------------------------------------------------------------------ #
fig = plt.figure(figsize=(7, 5))
ax  = fig.add_subplot(111, projection="3d")

def draw_frame(k):
    ax.clear()
    Z = frames[k]   # shape (N, N), Z[j,i] = u(x_i, y_j)
    surf = ax.plot_surface(
        X, Y, Z,
        cmap="RdBu_r",
        vmin=Z_MIN, vmax=Z_MAX,
        linewidth=0, antialiased=False,
        rcount=80, ccount=80        # resolution of the surface plot
    )
    ax.set_zlim(Z_MIN, Z_MAX)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("u [m]")
    ax.set_title(f"Sheet vibration — frame {k}/{nframes}")
    ax.view_init(elev=30, azim=-60)
    return surf,

# ------------------------------------------------------------------ #
#  Render animation                                                    #
# ------------------------------------------------------------------ #
print("Rendering frames...")
ani = animation.FuncAnimation(
    fig, draw_frame,
    frames=nframes,
    interval=1000 // 30,
    blit=False
)

print(f"Saving {OUTPUT_GIF} ...")
ani.save(OUTPUT_GIF, writer="pillow", fps=10)
print("Done.")
