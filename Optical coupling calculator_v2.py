import matplotlib
matplotlib.use('TkAgg')

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Ellipse

# Import your Tab 3 module
from fresnel_stack import FresnelStackTab

# --- Core Physics Functions ---

def calculate_fresnel_transmission(n_a, n_b):
    if n_a == 0 or n_b == 0: return 0
    R = ((n_a - n_b) / (n_a + n_b)) ** 2
    return 1 - R

def eta_M(Dx1, Dx2, Dy1, Dy2):
    if Dx1 == 0 or Dx2 == 0 or Dy1 == 0 or Dy2 == 0: return 0
    term_x = 2 / ((Dx1 / Dx2) + (Dx2 / Dx1))
    term_y = 2 / ((Dy1 / Dy2) + (Dy2 / Dy1))
    return term_x * term_y

def eta_delta_x(delta_x, Dx1, Dx2):
    if Dx1 == 0 or Dx2 == 0: return 0
    x_e = np.sqrt((Dx1**2 + Dx2**2) / 8)
    return np.exp(-(delta_x / x_e) ** 2)

def eta_delta_y(delta_y, Dy1, Dy2):
    if Dy1 == 0 or Dy2 == 0: return 0
    y_e = np.sqrt((Dy1**2 + Dy2**2) / 8)
    return np.exp(-(delta_y / y_e) ** 2)

def eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength):
    if Dx1 == 0 or Dx2 == 0 or Dy1 == 0 or Dy2 == 0 or wavelength == 0: return 0
    k = 2 * np.pi / wavelength
    # Rayleigh ranges (approximate for asymmetric beams)
    zRx = (np.pi * Dx1 * Dx2) / wavelength # Approximation logic from paper context usually
    
    # Using the specific overlap formula from the initial script:
    term_x = (Dx1 / Dx2 + Dx2 / Dx1) ** 2 + ((8 / (k * Dx1 * Dx2)) ** 2) * (delta_z ** 2)
    term_x = 2 / np.sqrt(term_x)
    
    term_y = (Dy1 / Dy2 + Dy2 / Dy1) ** 2 + ((8 / (k * Dy1 * Dy2)) ** 2) * (delta_z ** 2)
    term_y = 2 / np.sqrt(term_y)
    
    return term_x * term_y

def to_dB(value):
    if value <= 1e-9: return -99.99 
    return 10 * np.log10(value)

# --- Main App Logic ---

def calculate_results(*args):
    """Performs calculation AND updates plots."""
    try:
        # 1. Get Inputs
        Dx1 = float(entry_Dx1_var.get())
        Dx2 = float(entry_Dx2_var.get())
        Dy1 = float(entry_Dy1_var.get())
        Dy2 = float(entry_Dy2_var.get())
        delta_x = float(entry_delta_x_var.get())
        delta_y = float(entry_delta_y_var.get())
        delta_z = float(entry_delta_z_var.get())
        wavelength = float(entry_wavelength_var.get())
        n1 = float(entry_n1.get())
        n2 = float(entry_n2.get())
        
        z_sweep_max = float(entry_z_sweep_max_var.get())

        # 2. Calculations
        eta_m = eta_M(Dx1, Dx2, Dy1, Dy2)
        eta_dx = eta_delta_x(delta_x, Dx1, Dx2)
        eta_dy = eta_delta_y(delta_y, Dy1, Dy2)
        eta_zm = eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength)

        # Fresnel
        eta_fresnel_val = 1.0
        if var_add_fresnel.get():
            if delta_z == 0:
                eta_fresnel_val = calculate_fresnel_transmission(n1, n2)
            else:
                eta_fresnel_val = calculate_fresnel_transmission(n1, 1.0) * calculate_fresnel_transmission(1.0, n2)

        eta_total = eta_dx * eta_dy * eta_zm * eta_fresnel_val

        # 3. Update Text Results
        if var_dB.get():
            result_eta_M.set(f"{to_dB(eta_m):.4f}")
            result_eta_delta_x.set(f"{to_dB(eta_dx):.4f}")
            result_eta_delta_y.set(f"{to_dB(eta_dy):.4f}")
            result_eta_ZM.set(f"{to_dB(eta_zm):.4f}")
            result_eta_fresnel.set(f"{to_dB(eta_fresnel_val):.4f}" if var_add_fresnel.get() else "Off")
            result_eta_total.set(f"{to_dB(eta_total):.4f}")
        else:
            result_eta_M.set(f"{eta_m:.4f}")
            result_eta_delta_x.set(f"{eta_dx:.4f}")
            result_eta_delta_y.set(f"{eta_dy:.4f}")
            result_eta_ZM.set(f"{eta_zm:.4f}")
            result_eta_fresnel.set(f"{eta_fresnel_val:.4f}" if var_add_fresnel.get() else "Off")
            result_eta_total.set(f"{eta_total:.4f}")

        # 4. Trigger Plot Updates
        update_cross_section_plot(Dx1, Dy1, Dx2, Dy2, delta_x, delta_y)
        update_2d_heatmap(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength, delta_x, delta_y)
        update_1d_il_plot(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength, delta_x, delta_y)
        update_z_sweep_plot(Dx1, Dx2, Dy1, Dy2, wavelength, z_sweep_max, delta_z)

    except ValueError:
        pass

# --- Plotting Functions ---

def update_cross_section_plot(Dx1, Dy1, Dx2, Dy2, dx, dy):
    ax_cs.clear()
    
    # Port 1 (Reference at 0,0) - Blue
    e1 = Ellipse((0, 0), width=Dx1, height=Dy1, angle=0, 
                 edgecolor='blue', facecolor='blue', alpha=0.3, label='Port 1')
    e1_outline = Ellipse((0, 0), width=Dx1, height=Dy1, angle=0, 
                 edgecolor='blue', facecolor='none', lw=2)
    
    # Port 2 (Misaligned at dx, dy) - Red
    e2 = Ellipse((dx, dy), width=Dx2, height=Dy2, angle=0, 
                 edgecolor='red', facecolor='red', alpha=0.3, label='Port 2')
    e2_outline = Ellipse((dx, dy), width=Dx2, height=Dy2, angle=0, 
                 edgecolor='red', facecolor='none', lw=2, linestyle='--')
    
    ax_cs.add_patch(e1)
    ax_cs.add_patch(e1_outline)
    ax_cs.add_patch(e2)
    ax_cs.add_patch(e2_outline)
    
    max_dim = max(Dx1, Dx2, Dy1, Dy2, abs(dx)*2, abs(dy)*2) * 1.5
    if max_dim == 0: max_dim = 10
    
    ax_cs.set_xlim(-max_dim/2 + dx/2, max_dim/2 + dx/2)
    ax_cs.set_ylim(-max_dim/2 + dy/2, max_dim/2 + dy/2)
    ax_cs.set_aspect('equal')
    ax_cs.grid(True, linestyle=':', alpha=0.6)
    ax_cs.set_xlabel(label_CS_xlabel.get(), fontsize=7)
    ax_cs.set_ylabel(label_CS_ylabel.get(), fontsize=7)
    ax_cs.set_title(label_CS_title.get(), fontsize=9, pad=3)
    canvas_cs.draw()

def update_2d_heatmap(Dx1, Dx2, Dy1, Dy2, dz, wl, input_dx, input_dy):
    # 1. Clear the ENTIRE figure, not just the axis.
    # This automatically removes the old colorbar and axis, avoiding the crash.
    fig_hm.clear()
    
    # 2. Re-add the subplot fresh
    ax = fig_hm.add_subplot(111)
    
    # --- CALCULATION LOGIC ---
    max_range = max(abs(input_dx), abs(input_dy)) * 2
    if max_range < 5: max_range = 10 
    
    X = np.linspace(-max_range, max_range, 50)
    Y = np.linspace(-max_range, max_range, 50)
    X_g, Y_g = np.meshgrid(X, Y)
    Z_g = np.zeros_like(X_g)
    
    const_ZM = eta_ZM(Dx1, Dx2, Dy1, Dy2, dz, wl)
    
    for i in range(X_g.shape[0]):
        for j in range(X_g.shape[1]):
            val = eta_delta_x(X_g[i,j], Dx1, Dx2) * eta_delta_y(Y_g[i,j], Dy1, Dy2) * const_ZM
            Z_g[i,j] = to_dB(val)
            
    # --- PLOTTING ---
    contour = ax.contourf(X_g, Y_g, Z_g, levels=np.linspace(-10, 0, 50), cmap='magma')
    
    # Add colorbar (No need to check for existence or remove, the figure is fresh)
    cbar = fig_hm.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label_HM_ylabel.get(), fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    
    # Contour lines
    c1 = ax.contour(X_g, Y_g, Z_g, levels=[-1], colors='white', linewidths=0.8, linestyles='--')
    ax.clabel(c1, inline=True, fmt='-1 dB', fontsize=6, colors='white', manual=False)

    c2 = ax.contour(X_g, Y_g, Z_g, levels=[-3], colors='white', linewidths=0.8, linestyles=':')
    ax.clabel(c2, inline=True, fmt='-3 dB', fontsize=6, colors='white', manual=False)

    ax.scatter([input_dx], [input_dy], color='cyan', marker='x', s=50, label='Current')
    
    # Labels
    ax.set_xlabel(label_HM_xlabel.get(), fontsize=7)
    ax.set_ylabel(label_HM_ylabel.get(), fontsize=7)
    ax.set_title(label_HM_title.get(), fontsize=9, pad=3)
    
    # 3. Draw
    canvas_hm.draw()

def update_1d_il_plot(Dx1, Dx2, Dy1, Dy2, dz, wl, input_dx, input_dy):
    ax_il.clear()
    
    max_delta = max(abs(input_dx), abs(input_dy)) * 1.5
    if max_delta < 5: max_delta = 10
    
    delta_vals = np.linspace(-max_delta, max_delta, 100)
    
    # Pre-calc constant ZM
    c_ZM = eta_ZM(Dx1, Dx2, Dy1, Dy2, dz, wl)
    
    # X Sweep (y=0)
    il_x = [to_dB(eta_delta_x(d, Dx1, Dx2) * eta_delta_y(0, Dy1, Dy2) * c_ZM) for d in delta_vals]
    
    # Y Sweep (x=0)
    il_y = [to_dB(eta_delta_x(0, Dx1, Dx2) * eta_delta_y(d, Dy1, Dy2) * c_ZM) for d in delta_vals]
    
    ax_il.plot(delta_vals, il_x, label=label_1D_dx_label.get(), color='blue')
    ax_il.plot(delta_vals, il_y, label=label_1D_dy_label.get(), color='green', linestyle='--')
    
    # Current Points
    curr_x = to_dB(eta_delta_x(input_dx, Dx1, Dx2) * eta_delta_y(0, Dy1, Dy2) * c_ZM)
    curr_y = to_dB(eta_delta_x(0, Dx1, Dx2) * eta_delta_y(input_dy, Dy1, Dy2) * c_ZM)
    
    ax_il.scatter([input_dx], [curr_x], color='blue', zorder=5)
    ax_il.scatter([input_dy], [curr_y], color='green', zorder=5)

    ax_z.set_xlabel(label_Z_xlabel.get(), fontsize=7)
    ax_z.set_ylabel(label_Z_ylabel.get(), fontsize=7)
    ax_z.set_title(label_Z_title.get(), fontsize=9, pad=3)
    ax_il.grid(True, alpha=0.5)
    ax_il.legend(fontsize=6)
    canvas_il.draw()

def update_z_sweep_plot(Dx1, Dx2, Dy1, Dy2, wl, z_max, current_z):
    ax_z.clear()
    
    if z_max <= 0: z_max = 50
    z_vals = np.linspace(0, z_max, 100)
    
    # Assume lateral misalignment is 0 for the Z sweep curve to show ideal Z behavior,
    # OR we can include current lateral misalignment. 
    # Typically Z-curves show the BEST CASE for Z, so we use lateral=0.
    # However, to be consistent with "Insertion Loss", we should probably include the fixed lateral loss.
    # Let's include the fixed lateral loss so the curve passes through the actual operating point.
    
    lateral_loss_factor = eta_delta_x(0, Dx1, Dx2) * eta_delta_y(0, Dy1, Dy2) # Ideal lateral
    # NOTE: If we want to show the specific loss at current (dx, dy) vs Z:
    # lateral_loss_factor = eta_delta_x(entry_delta_x_var.get()...)
    # Standard practice is usually ideal lateral, but let's stick to the misalignment provided to be strict.
    
    il_z = []
    for z in z_vals:
        val = lateral_loss_factor * eta_ZM(Dx1, Dx2, Dy1, Dy2, z, wl)
        il_z.append(to_dB(val))
        
    ax_z.plot(z_vals, il_z, color='purple', label='IL vs Δz')
    
    # Current Z point
    curr_val = lateral_loss_factor * eta_ZM(Dx1, Dx2, Dy1, Dy2, current_z, wl)
    ax_z.scatter([current_z], [to_dB(curr_val)], color='purple', zorder=5)
    
    ax_z.set_xlabel("Δz (μm)", fontsize=7)
    ax_z.set_ylabel("IL (dB)", fontsize=7)
    ax_z.set_title("Longitudinal Loss", fontsize=9, pad=3)
    ax_z.grid(True, alpha=0.5)
    canvas_z.draw()

# --- UI Helper Functions ---

def toggle_symmetry_Dx():
    if var_symmetry_Dx.get():
        entry_Dy1.config(state='disabled')
        entry_Dy1_var.set(entry_Dx1_var.get())
    else:
        entry_Dy1.config(state='normal')
    calculate_results()

def toggle_symmetry_Dy():
    if var_symmetry_Dy.get():
        entry_Dy2.config(state='disabled')
        entry_Dy2_var.set(entry_Dx2_var.get())
    else:
        entry_Dy2.config(state='normal')
    calculate_results()

def sync_Dy1(*args):
    if var_symmetry_Dx.get(): entry_Dy1_var.set(entry_Dx1_var.get())
    calculate_results()

def sync_Dy2(*args):
    if var_symmetry_Dy.get(): entry_Dy2_var.set(entry_Dx2_var.get())
    calculate_results()

def open_label_settings():
    label_window = tk.Toplevel(root)
    label_window.title("Edit Plot Labels")
    label_window.geometry("350x300")
    
    # Create Tabs inside the popup
    nb = ttk.Notebook(label_window)
    nb.pack(fill='both', expand=True, padx=5, pady=5)
    
    tab_1d = ttk.Frame(nb, padding=10)
    tab_cs = ttk.Frame(nb, padding=10)
    tab_hm = ttk.Frame(nb, padding=10)
    tab_z  = ttk.Frame(nb, padding=10)
    
    nb.add(tab_1d, text='1D Plot')
    nb.add(tab_cs, text='Cross Sect')
    nb.add(tab_hm, text='Heatmap')
    nb.add(tab_z,  text='Z Sweep')
    
    # Helper to build rows
    def add_rows(parent, items):
        for i, (text, var) in enumerate(items):
            ttk.Label(parent, text=text).grid(row=i, column=0, sticky='w', pady=5)
            e = ttk.Entry(parent, textvariable=var, width=20)
            e.grid(row=i, column=1, sticky='w', padx=10, pady=5)
            e.bind("<KeyRelease>", lambda e: calculate_results())

    # Define items for each tab
    items_1d = [
        ("Title:", label_1D_title),
        ("X Axis:", label_1D_xlabel),
        ("Y Axis:", label_1D_ylabel),
        ("Legend X:", label_1D_dx_label),
        ("Legend Y:", label_1D_dy_label)
    ]
    
    items_cs = [
        ("Title:", label_CS_title),
        ("X Axis:", label_CS_xlabel),
        ("Y Axis:", label_CS_ylabel)
    ]
    
    items_hm = [
        ("Title:", label_HM_title),
        ("X Axis:", label_HM_xlabel),
        ("Y Axis:", label_HM_ylabel)
    ]
    
    items_z = [
        ("Title:", label_Z_title),
        ("X Axis:", label_Z_xlabel),
        ("Y Axis:", label_Z_ylabel)
    ]

    add_rows(tab_1d, items_1d)
    add_rows(tab_cs, items_cs)
    add_rows(tab_hm, items_hm)
    add_rows(tab_z, items_z)

# --- Main Setup ---

root = tk.Tk()
root.title("Optical Coupling Calculator")
root.geometry("1100x700") # Wider for 4 plots

notebook = ttk.Notebook(root)
notebook.pack(fill='both', expand=True)

tab1 = ttk.Frame(notebook)
tab2 = ttk.Frame(notebook)
tab3 = FresnelStackTab(notebook)

notebook.add(tab1, text='Coupling Calculator')
notebook.add(tab2, text='Converter')
notebook.add(tab3, text='Fresnel Stack')

# ================= TAB 1 LAYOUT =================

# Split Tab 1: Left (Inputs) | Right (Plots)
frame_left = ttk.Frame(tab1, padding=10)
frame_left.grid(row=0, column=0, sticky="nsew")

frame_right = ttk.Frame(tab1, padding=10)
frame_right.grid(row=0, column=1, sticky="nsew")

tab1.columnconfigure(1, weight=1)
tab1.rowconfigure(0, weight=1)

# --- INPUT VARIABLES ---
entry_Dx1_var = tk.StringVar(value="10")
entry_Dx2_var = tk.StringVar(value="10")
entry_Dy1_var = tk.StringVar(value="10")
entry_Dy2_var = tk.StringVar(value="10")
entry_delta_x_var = tk.StringVar(value="2")
entry_delta_y_var = tk.StringVar(value="0")
entry_delta_z_var = tk.StringVar(value="0")
entry_wavelength_var = tk.StringVar(value="1.55")
entry_z_sweep_max_var = tk.StringVar(value="50") # New Input

# Traces
for var in [entry_Dx1_var, entry_Dx2_var, entry_Dy1_var, entry_Dy2_var, 
            entry_delta_x_var, entry_delta_y_var, entry_delta_z_var, 
            entry_wavelength_var, entry_z_sweep_max_var]:
    var.trace_add('write', calculate_results)

# --- LEFT COLUMN ---
# Geometry
ttk.Label(frame_left, text="Dx1 (μm):").grid(row=0, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_Dx1_var, width=8).grid(row=0, column=1)
var_symmetry_Dx = tk.BooleanVar()
ttk.Checkbutton(frame_left, text="Sym", variable=var_symmetry_Dx, command=toggle_symmetry_Dx).grid(row=0, column=2)

ttk.Label(frame_left, text="Dx2 (μm):").grid(row=1, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_Dx2_var, width=8).grid(row=1, column=1)
var_symmetry_Dy = tk.BooleanVar()
ttk.Checkbutton(frame_left, text="Sym", variable=var_symmetry_Dy, command=toggle_symmetry_Dy).grid(row=1, column=2)

ttk.Label(frame_left, text="Dy1 (μm):").grid(row=2, column=0, sticky='w')
entry_Dy1 = ttk.Entry(frame_left, textvariable=entry_Dy1_var, width=8)
entry_Dy1.grid(row=2, column=1)

ttk.Label(frame_left, text="Dy2 (μm):").grid(row=3, column=0, sticky='w')
entry_Dy2 = ttk.Entry(frame_left, textvariable=entry_Dy2_var, width=8)
entry_Dy2.grid(row=3, column=1)

ttk.Separator(frame_left, orient='horizontal').grid(row=4, column=0, columnspan=3, sticky='ew', pady=5)

ttk.Label(frame_left, text="δx (μm):").grid(row=5, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_delta_x_var, width=8).grid(row=5, column=1)

ttk.Label(frame_left, text="δy (μm):").grid(row=6, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_delta_y_var, width=8).grid(row=6, column=1)

ttk.Label(frame_left, text="δz (μm):").grid(row=7, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_delta_z_var, width=8).grid(row=7, column=1)

ttk.Label(frame_left, text="λ (μm):").grid(row=8, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_wavelength_var, width=8).grid(row=8, column=1)

# New Z Sweep Input
ttk.Label(frame_left, text="Z Sweep Max:").grid(row=9, column=0, sticky='w')
ttk.Entry(frame_left, textvariable=entry_z_sweep_max_var, width=8).grid(row=9, column=1)
ttk.Label(frame_left, text="μm").grid(row=9, column=2, sticky='w')


ttk.Separator(frame_left, orient='horizontal').grid(row=10, column=0, columnspan=3, sticky='ew', pady=5)

# Fresnel Inputs
ttk.Label(frame_left, text="n1:").grid(row=11, column=0, sticky='w')
entry_n1 = ttk.Entry(frame_left, width=8); entry_n1.insert(0, "3.2"); entry_n1.grid(row=11, column=1)
ttk.Label(frame_left, text="n2:").grid(row=12, column=0, sticky='w')
entry_n2 = ttk.Entry(frame_left, width=8); entry_n2.insert(0, "1.53"); entry_n2.grid(row=12, column=1)

# Results
start_res = 13
result_eta_M = tk.StringVar()
result_eta_delta_x = tk.StringVar()
result_eta_delta_y = tk.StringVar()
result_eta_ZM = tk.StringVar()
result_eta_fresnel = tk.StringVar()
result_eta_total = tk.StringVar()

labels_res = [("η_M:", result_eta_M), ("η_δx:", result_eta_delta_x), 
              ("η_δy:", result_eta_delta_y), ("η_ZM:", result_eta_ZM),
              ("η_Fres:", result_eta_fresnel), ("Total:", result_eta_total)]

for i, (txt, var) in enumerate(labels_res):
    ttk.Label(frame_left, text=txt, font=('Arial', 9, 'bold')).grid(row=start_res+i, column=0, sticky='w', pady=1)
    ttk.Label(frame_left, textvariable=var).grid(row=start_res+i, column=1, sticky='w', pady=1)

# Buttons
btn_row = start_res + len(labels_res) + 1
var_dB = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_left, text="dB", variable=var_dB, command=calculate_results).grid(row=btn_row, column=0)
var_add_fresnel = tk.BooleanVar()
ttk.Checkbutton(frame_left, text="+Fresnel", variable=var_add_fresnel, command=calculate_results).grid(row=btn_row, column=1, columnspan=2)

ttk.Button(frame_left, text="Labels", command=open_label_settings).grid(row=btn_row+1, column=0, columnspan=2, pady=5)

# --- RIGHT COLUMN (2x2 PLOTS) ---

frame_right.columnconfigure(0, weight=1)
frame_right.columnconfigure(1, weight=1)
frame_right.rowconfigure(0, weight=1)
frame_right.rowconfigure(1, weight=1)

# Variables for Labels
label_1D_title = tk.StringVar(value="1D Insertion Loss")
label_1D_xlabel = tk.StringVar(value="Misalignment (μm)")
label_1D_ylabel = tk.StringVar(value="IL (dB)")
label_1D_dx_label = tk.StringVar(value="X Sweep")
label_1D_dy_label = tk.StringVar(value="Y Sweep")

# Cross-Section (Top-Left)
label_CS_title = tk.StringVar(value="Cross-Section")
label_CS_xlabel = tk.StringVar(value="x (μm)")
label_CS_ylabel = tk.StringVar(value="y (μm)")

# Heatmap (Top-Right)
label_HM_title = tk.StringVar(value="2D Heatmap (dB)")
label_HM_xlabel = tk.StringVar(value="δx (μm)")
label_HM_ylabel = tk.StringVar(value="δy (μm)")

# Z Sweep (Bottom-Right)
label_Z_title = tk.StringVar(value="Longitudinal Loss")
label_Z_xlabel = tk.StringVar(value="Δz (μm)")
label_Z_ylabel = tk.StringVar(value="IL (dB)")


# 1. Top-Left: Cross Section
fig_cs = plt.Figure(figsize=(3, 3), dpi=100)
ax_cs = fig_cs.add_subplot(111)
canvas_cs = FigureCanvasTkAgg(fig_cs, master=frame_right)
canvas_cs.get_tk_widget().grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

# 2. Top-Right: 2D Heatmap
fig_hm = plt.Figure(figsize=(3, 3), dpi=100)
ax_hm = fig_hm.add_subplot(111)
canvas_hm = FigureCanvasTkAgg(fig_hm, master=frame_right)
canvas_hm.get_tk_widget().grid(row=0, column=1, sticky="nsew", padx=2, pady=2)

# 3. Bottom-Left: 1D IL
fig_il = plt.Figure(figsize=(3, 3), dpi=100)
ax_il = fig_il.add_subplot(111)
fig_il.subplots_adjust(bottom=0.2, left=0.2)
canvas_il = FigureCanvasTkAgg(fig_il, master=frame_right)
canvas_il.get_tk_widget().grid(row=1, column=0, sticky="nsew", padx=2, pady=2)

# 4. Bottom-Right: Z Sweep
fig_z = plt.Figure(figsize=(3, 3), dpi=100)
ax_z = fig_z.add_subplot(111)
fig_z.subplots_adjust(bottom=0.2, left=0.2)
canvas_z = FigureCanvasTkAgg(fig_z, master=frame_right)
canvas_z.get_tk_widget().grid(row=1, column=1, sticky="nsew", padx=2, pady=2)

# Initialize
calculate_results()

# ================= TAB 2: CONVERTER =================
# (Placeholder for existing Tab 2 code)

# Percentage to dB Section
ttk.Label(tab2, text="Percentage to dB", font=('Arial', 11, 'bold')).grid(row=0, column=0, columnspan=2, sticky='w', padx=5, pady=(10,5))
ttk.Label(tab2, text="Percentage (%):").grid(row=1, column=0, sticky='w', padx=5, pady=2)
entry_percent_to_db = ttk.Entry(tab2)
entry_percent_to_db.insert(0, "100")
entry_percent_to_db.grid(row=1, column=1, padx=5, pady=2)

result_percent_to_db = tk.StringVar()
ttk.Label(tab2, textvariable=result_percent_to_db).grid(row=2, column=0, columnspan=2, padx=5, pady=2)

def calc_percent_to_db():
    try:
        percent = float(entry_percent_to_db.get())
        if percent <= 0:
            result_percent_to_db.set("Percentage must be > 0")
        else:
            db_value = 10 * np.log10(percent / 100)
            result_percent_to_db.set(f"Result: {db_value:.4f} dB")
    except ValueError:
        result_percent_to_db.set("Invalid input")

ttk.Button(tab2, text="Calculate", command=calc_percent_to_db).grid(row=1, column=2, padx=5, pady=2)

# Separator
ttk.Separator(tab2, orient='horizontal').grid(row=3, column=0, columnspan=3, sticky='ew', pady=10)

# dB to Percentage Section
ttk.Label(tab2, text="dB to Percentage", font=('Arial', 11, 'bold')).grid(row=4, column=0, columnspan=2, sticky='w', padx=5, pady=5)
ttk.Label(tab2, text="dB:").grid(row=5, column=0, sticky='w', padx=5, pady=2)
entry_db_to_percent = ttk.Entry(tab2)
entry_db_to_percent.insert(0, "0")
entry_db_to_percent.grid(row=5, column=1, padx=5, pady=2)

result_db_to_percent = tk.StringVar()
ttk.Label(tab2, textvariable=result_db_to_percent).grid(row=6, column=0, columnspan=2, padx=5, pady=2)

def calc_db_to_percent():
    try:
        db_value = float(entry_db_to_percent.get())
        percent = 100 * (10 ** (db_value / 10))
        result_db_to_percent.set(f"Result: {percent:.4f} %")
    except ValueError:
        result_db_to_percent.set("Invalid input")

ttk.Button(tab2, text="Calculate", command=calc_db_to_percent).grid(row=5, column=2, padx=5, pady=2)

# Separator
ttk.Separator(tab2, orient='horizontal').grid(row=7, column=0, columnspan=3, sticky='ew', pady=10)

# Current to Optical Power Section
ttk.Label(tab2, text="Current to Optical Power", font=('Arial', 11, 'bold')).grid(row=8, column=0, columnspan=2, sticky='w', padx=5, pady=5)
ttk.Label(tab2, text="Current (mA):").grid(row=9, column=0, sticky='w', padx=5, pady=2)
entry_current = ttk.Entry(tab2)
entry_current.insert(0, "0.001")
entry_current.grid(row=9, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Responsivity (A/W):").grid(row=10, column=0, sticky='w', padx=5, pady=2)
entry_responsivity = ttk.Entry(tab2)
entry_responsivity.insert(0, "0.8")
entry_responsivity.grid(row=10, column=1, padx=5, pady=2)

result_optical_power = tk.StringVar()
ttk.Label(tab2, textvariable=result_optical_power).grid(row=11, column=0, columnspan=2, padx=5, pady=2)

def calc_optical_power():
    try:
        current = float(entry_current.get())
        responsivity = float(entry_responsivity.get())
        if responsivity == 0:
            result_optical_power.set("Responsivity cannot be 0")
        else:
            power = current / responsivity
            power_dbm = 10 * np.log10(power)  # Convert to dBm
            result_optical_power.set(f"Result: {power_dbm:.2f} dBm, ({power:.4f} mW)")
    except ValueError:
        result_optical_power.set("Invalid input")

ttk.Button(tab2, text="Calculate", command=calc_optical_power).grid(row=9, column=2, rowspan=2, padx=5, pady=2)

# Separator
ttk.Separator(tab2, orient='horizontal').grid(row=12, column=0, columnspan=3, sticky='ew', pady=10)

# Half Angle Divergence Section
ttk.Label(tab2, text="Half Angle Divergence", font=('Arial', 11, 'bold')).grid(row=13, column=0, columnspan=2, sticky='w', padx=5, pady=5)
ttk.Label(tab2, text="Axial Distance, z (µm):").grid(row=14, column=0, sticky='w', padx=5, pady=2)
entry_z_dist = ttk.Entry(tab2)
entry_z_dist.insert(0, "10")
entry_z_dist.grid(row=14, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Beam Waist, ω₀ (µm):").grid(row=15, column=0, sticky='w', padx=5, pady=2)
entry_beam_waist = ttk.Entry(tab2)
entry_beam_waist.insert(0, "1")
entry_beam_waist.grid(row=15, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Wavelength, λ (µm):").grid(row=16, column=0, sticky='w', padx=5, pady=2)
entry_wl_div = ttk.Entry(tab2)
entry_wl_div.insert(0, "0.633")
entry_wl_div.grid(row=16, column=1, padx=5, pady=2)

result_divergence = tk.StringVar()
result_beam_diameter = tk.StringVar()
result_rayleigh = tk.StringVar()
result_rayleigh_half = tk.StringVar()
result_half_angle = tk.StringVar()

ttk.Label(tab2, text="Half Beam Diameter, ω(z):").grid(row=17, column=0, sticky='w', padx=5, pady=2)
ttk.Label(tab2, textvariable=result_beam_diameter).grid(row=17, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Radius of Curvature, R(z):").grid(row=18, column=0, sticky='w', padx=5, pady=2)
ttk.Label(tab2, textvariable=result_divergence).grid(row=18, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Rayleigh Range, Z_R:").grid(row=19, column=0, sticky='w', padx=5, pady=2)
ttk.Label(tab2, textvariable=result_rayleigh).grid(row=19, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Rayleigh Half Diameter, ω_R:").grid(row=20, column=0, sticky='w', padx=5, pady=2)
ttk.Label(tab2, textvariable=result_rayleigh_half).grid(row=20, column=1, padx=5, pady=2)

ttk.Label(tab2, text="Half Angle Divergence, θ (rad):").grid(row=21, column=0, sticky='w', padx=5, pady=2)
ttk.Label(tab2, textvariable=result_half_angle).grid(row=21, column=1, padx=5, pady=2)

def calc_divergence():
    try:
        z = float(entry_z_dist.get()) *1e-3  # mm
        w0 = float(entry_beam_waist.get()) *1e-3 # mm
        wavelength = float(entry_wl_div.get())  # µm
        
        # Convert wavelength to mm for consistency
        wavelength_mm = wavelength / 1000
        
        # Calculate Rayleigh range
        z_R = (np.pi * w0**2) / wavelength_mm
        
        # Calculate beam diameter at distance z
        w_z = w0 * np.sqrt(1 + (z / z_R)**2)
        
        # Calculate radius of curvature
        if z != 0:
            R_z = z * (1 + (z_R / z)**2)
        else:
            R_z = float('inf')
        
        # Calculate Rayleigh half diameter
        w_R = np.sqrt(2) * w0
        
        # Calculate half angle divergence in radians and degrees
        theta_rad = wavelength_mm / (np.pi * w0)
        theta_deg = np.degrees(theta_rad)
        theta_mrad = theta_rad * 1000
        
        result_beam_diameter.set(f"{w_z:.4f} mm")
        if R_z == float('inf'):
            result_divergence.set("∞ mm")
        else:
            result_divergence.set(f"{R_z:.4f} mm")
        result_rayleigh.set(f"{z_R:.4f} mm")
        result_rayleigh_half.set(f"{w_R:.4f} mm")
        result_half_angle.set(f"{theta_rad:.6f} rad, ({theta_deg:.4f}°)")
        
    except ValueError:
        result_beam_diameter.set("Invalid input")
        result_divergence.set("Invalid input")
        result_rayleigh.set("Invalid input")
        result_rayleigh_half.set("Invalid input")
        result_half_angle.set("Invalid input")

ttk.Button(tab2, text="Calculate", command=calc_divergence).grid(row=14, column=2, rowspan=3, padx=5, pady=2)

# Separator
ttk.Separator(tab2, orient='horizontal').grid(row=22, column=0, columnspan=3, sticky='ew', pady=10)

# dBm / mW Conversion Section
ttk.Label(tab2, text="dBm / mW Conversion", font=('Arial', 11, 'bold')).grid(
    row=23, column=0, columnspan=2, sticky='w', padx=5, pady=5
)

# dBm → mW
ttk.Label(tab2, text="Power (dBm):").grid(row=24, column=0, sticky='w', padx=5, pady=2)
entry_dbm = ttk.Entry(tab2)
entry_dbm.insert(0, "0")
entry_dbm.grid(row=24, column=1, padx=5, pady=2)

result_mw = tk.StringVar()
ttk.Label(tab2, textvariable=result_mw).grid(row=25, column=0, columnspan=2, padx=5, pady=2)

def conv_dbm_to_mw():
    try:
        dbm = float(entry_dbm.get())
        mw = 10**(dbm / 10)
        result_mw.set(f"{mw:.4f} mW")
    except ValueError:
        result_mw.set("Invalid input")

ttk.Button(tab2, text="Convert", command=conv_dbm_to_mw).grid(
    row=24, column=2, padx=5, pady=2
)

# mW → dBm
ttk.Label(tab2, text="Power (mW):").grid(row=26, column=0, sticky='w', padx=5, pady=2)
entry_mw = ttk.Entry(tab2)
entry_mw.insert(0, "1")
entry_mw.grid(row=26, column=1, padx=5, pady=2)

result_dbm = tk.StringVar()
ttt = ttk.Label(tab2, textvariable=result_dbm).grid(row=27, column=0, columnspan=2, padx=5, pady=2)

def conv_mw_to_dbm():
    try:
        mw = float(entry_mw.get())
        if mw <= 0:
            result_dbm.set("mW must be > 0")
            return
        dbm = 10 * np.log10(mw)
        result_dbm.set(f"{dbm:.2f} dBm")
    except ValueError:
        result_dbm.set("Invalid input")

ttk.Button(tab2, text="Convert", command=conv_mw_to_dbm).grid(
    row=26, column=2, padx=5, pady=2
)


# --- Cleanup ---
def on_closing():
    plt.close('all')
    root.quit()
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)

root.mainloop()
