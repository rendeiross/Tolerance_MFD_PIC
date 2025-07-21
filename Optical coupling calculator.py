import matplotlib
matplotlib.use('TkAgg')

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import math as math

#The equations used were retrieved from an open source paper: https://doi.org/10.1364/AO.540682


def eta_M(Dx1, Dx2, Dy1, Dy2):
    """Calculate η_M based on the given diameters."""
    if Dx1 == 0 or Dx2 == 0 or Dy1 == 0 or Dy2 == 0:
        return 0
    term_x = 2 / ((Dx1 / Dx2) + (Dx2 / Dx1))
    term_y = 2 / ((Dy1 / Dy2) + (Dy2 / Dy1))
    return term_x * term_y

def eta_delta_x(delta_x, Dx1, Dx2):
    """Calculate η_δx using the misalignment parameter."""
    if Dx1 == 0 or Dx2 == 0:
        return 0
    x_e = np.sqrt((Dx1**2 + Dx2**2) / 8)
    return np.exp(-(delta_x / x_e) ** 2)

def eta_delta_y(delta_y, Dy1, Dy2):
    """Calculate η_δy using the misalignment parameter."""
    if Dy1 == 0 or Dy2 == 0:
        return 0
    y_e = np.sqrt((Dy1**2 + Dy2**2) / 8)
    return np.exp(-(delta_y / y_e) ** 2)

def eta_delta_teta(delta_teta, Dx1, Dx2, Dy1, Dy2, wavelength,n):
    """Calculate η_δteta using the misalignment parameter."""
    D1 = np.sqrt(Dx1**2+Dy1**2)
    D2 = np.sqrt(Dx2**2+Dy2**2)
    teta = (2*np.sqrt(2)*wavelength/np.pi*n)*np.sqrt(1/(D1**2+D2**2))
    efficiency = np.exp((-(np.tan(delta_teta)/teta)**2))
    return efficiency

def eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength):
    """Calculate η_ZM based on the given diameters and misalignment parameter delta_z."""
    if Dx1 == 0 or Dx2 == 0 or Dy1 == 0 or Dy2 == 0 or wavelength == 0:
        return 0
    k = 2 * np.pi / wavelength  # Wavenumber

    term_x = (Dx1 / Dx2 + Dx2 / Dx1) ** 2
    term_x += ((8 / (k * Dx1 * Dx2)) ** 2) * (delta_z ** 2)
    term_x = 2 / np.sqrt(term_x)

    term_y = (Dy1 / Dy2 + Dy2 / Dy1) ** 2
    term_y += ((8 / (k * Dy1 * Dy2)) ** 2) * (delta_z ** 2)
    term_y = 2 / np.sqrt(term_y)

    return term_x * term_y

def to_dB(value):
    """Convert efficiency to dB."""
    with np.errstate(divide='ignore'):
        return 10 * np.log10(value)

def calculate_results():
    try:
        Dx1 = float(entry_Dx1.get())
        Dx2 = float(entry_Dx2.get())
        Dy1 = float(entry_Dy1.get())
        Dy2 = float(entry_Dy2.get())
        delta_x = float(entry_delta_x.get())
        delta_teta = math.radians(float((entry_delta_teta.get())))
        delta_phi = math.radians(float(entry_delta_phi.get()))
        n = float(entry_n.get())
        delta_y = float(entry_delta_y.get())
        delta_z = float(entry_delta_z.get())
        wavelength = float(entry_wavelength.get())

        eta_m_value = round(eta_M(Dx1, Dx2, Dy1, Dy2), 4)
        eta_dx_value = round(eta_delta_x(delta_x, Dx1, Dx2), 4)
        eta_dy_value = round(eta_delta_y(delta_y, Dy1, Dy2), 4)
        eta_delta_teta_value = eta_delta_teta(delta_teta, Dx1, Dx2, Dy1, Dy2, wavelength,n)#round(eta_delta_teta(delta_teta, Dx1, Dx2, Dy1, Dy2, wavelength,n) 4)
        eta_delta_phi_value = eta_delta_teta(delta_phi, Dx1, Dx2, Dy1, Dy2, wavelength,n)#round(eta_delta_teta(delta_phi, Dx1, Dx2, Dy1, Dy2, wavelength,n) 4)
        eta_zm_value = round(eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength), 4)
        eta_total_value = round(eta_dx_value * eta_dy_value * eta_zm_value * eta_delta_teta_value * eta_delta_phi_value, 4)

        if var_dB.get():
            result_eta_M.set(round(to_dB(eta_m_value), 4))
            result_eta_delta_x.set(round(to_dB(eta_dx_value), 4))
            result_eta_delta_y.set(round(to_dB(eta_dy_value), 4))
            result_eta_delta_teta.set(round(to_dB(eta_delta_teta_value), 4))
            result_eta_delta_phi.set(round(to_dB(eta_delta_phi_value), 4))
            result_eta_ZM.set(round(to_dB(eta_zm_value), 4))
            result_eta_total.set(round(to_dB(eta_total_value), 4))
        else:
            result_eta_M.set(eta_m_value)
            result_eta_delta_x.set(eta_dx_value)
            result_eta_delta_y.set(eta_dy_value)
            result_eta_delta_teta.set(eta_delta_teta_value)
            result_eta_delta_phi.set(eta_delta_phi_value)
            result_eta_ZM.set(eta_zm_value)
            result_eta_total.set(eta_total_value)
    except ValueError:
        result_eta_M.set("Invalid input")
        result_eta_delta_x.set("Invalid input")
        result_eta_delta_y.set("Invalid input")
        result_eta_delta_teta.set("Invalid input")
        result_eta_delta_phi.set("Invalid input")
        result_eta_ZM.set("Invalid input")
        result_eta_total.set("Invalid input")

def toggle_symmetry_Dx():
    if var_symmetry_Dx.get():
        entry_Dy1.config(state='disabled')
        entry_Dy1_var.set(entry_Dx1_var.get())
    else:
        entry_Dy1.config(state='normal')

def toggle_symmetry_Dy():
    if var_symmetry_Dy.get():
        entry_Dy2.config(state='disabled')
        entry_Dy2_var.set(entry_Dx2_var.get())
    else:
        entry_Dy2.config(state='normal')

def sync_Dy1(*args):
    if var_symmetry_Dx.get():
        entry_Dy1_var.set(entry_Dx1_var.get())

def sync_Dy2(*args):
    if var_symmetry_Dy.get():
        entry_Dy2_var.set(entry_Dx2_var.get())


def update_results():
    calculate_results()

def generate_heatmap():
    try:
        Dx1 = float(entry_Dx1.get())
        Dx2 = float(entry_Dx2.get())
        Dy1 = float(entry_Dy1.get())
        Dy2 = float(entry_Dy2.get())
        delta_z = float(entry_delta_z.get())
        wavelength = float(entry_wavelength.get())
        
        dx_range = float(entry_delta_x.get())
        dy_range = float(entry_delta_y.get())
        
        X_d = np.linspace(-dx_range, dx_range, 100)
        Y_d = np.linspace(-dy_range, dy_range, 100)
        X_d, Y_d = np.meshgrid(X_d, Y_d)
        
        efficiency_losses = np.zeros_like(X_d)
        
        for i in range(X_d.shape[0]):
            for j in range(X_d.shape[1]):
                efficiency_losses[i,j] = eta_delta_x(X_d[i,j], Dx1, Dx2) * \
                                         eta_delta_y(Y_d[i,j], Dy1, Dy2) * \
                                         eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength)
        
        efficiency_losses_db = to_dB(efficiency_losses)
        
        plt.figure()
        contour = plt.contourf(X_d, Y_d, efficiency_losses_db, levels=np.linspace(-6, 0, 500), cmap='magma')
        plt.colorbar(contour, label='Insertion Loss (dB)')
        
        # Add -1 dB and -3 dB contour lines
        contour_line1 = plt.contour(X_d, Y_d, efficiency_losses_db, levels=[-1], colors='black', linewidths=1, linestyles='dashed')
        plt.clabel(contour_line1, fmt='-1 dB', inline=True, fontsize=8)
        
        
        contour_line1 = plt.contour(X_d, Y_d, efficiency_losses_db, levels=[-2], colors='black', linewidths=1, linestyles='dashed')
        plt.clabel(contour_line1, fmt='-2 dB', inline=True, fontsize=8)

        contour_line3 = plt.contour(X_d, Y_d, efficiency_losses_db, levels=[-3], colors='black', linewidths=1, linestyles='dotted')
        plt.clabel(contour_line3, fmt='-3 dB', inline=True, fontsize=8)
        
        plt.xlabel(label_2D_xlabel.get())
        plt.ylabel(label_2D_ylabel.get())
        plt.title(label_2D_title.get())
        
        plt.show()
        
    except ValueError:
        print("Invalid input for heatmap generation")

def generate_1D_IL_map():
    try:
        Dx1 = float(entry_Dx1.get())
        Dx2 = float(entry_Dx2.get())
        Dy1 = float(entry_Dy1.get())
        Dy2 = float(entry_Dy2.get())
        delta_z = float(entry_delta_z.get())
        wavelength = float(entry_wavelength.get())

        max_delta = max(float(entry_delta_x.get()), float(entry_delta_y.get()))
        delta_values = np.linspace(-max_delta, max_delta, 200)

        il_x = []
        il_y = []

        for delta in delta_values:
            eta_total_x =  eta_delta_x(delta, Dx1, Dx2) * eta_delta_y(0, Dy1, Dy2) * eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength)
            eta_total_y = eta_delta_x(0, Dx1, Dx2) * eta_delta_y(delta, Dy1, Dy2) * eta_ZM(Dx1, Dx2, Dy1, Dy2, delta_z, wavelength)
            il_x.append(to_dB(eta_total_x))
            il_y.append(to_dB(eta_total_y))

        plt.figure()
        plt.plot(delta_values, il_x, label=label_1D_dx_label.get(), color='blue')
        plt.plot(delta_values, il_y, label=label_1D_dy_label.get(), color='green')
        plt.axhline(-1, color='gray', linestyle='--', linewidth=0.8)
        plt.axhline(-3, color='gray', linestyle='--', linewidth=0.8)
        plt.xlabel(label_1D_xlabel.get())
        plt.ylabel(label_1D_ylabel.get())
        plt.title(label_1D_title.get())
        plt.legend()
        plt.grid(True)
        plt.show()
        
    except ValueError:
        print("Invalid input for 1D IL plot")

def open_label_settings():
    label_window = tk.Toplevel(root)
    label_window.title("Edit Plot Labels")

    # 2D plot labels
    ttk.Label(label_window, text="2D Plot Title").grid(row=0, column=0)
    ttk.Entry(label_window, textvariable=label_2D_title).grid(row=0, column=1)

    ttk.Label(label_window, text="2D X Label").grid(row=1, column=0)
    ttk.Entry(label_window, textvariable=label_2D_xlabel).grid(row=1, column=1)

    ttk.Label(label_window, text="2D Y Label").grid(row=2, column=0)
    ttk.Entry(label_window, textvariable=label_2D_ylabel).grid(row=2, column=1)

    # 1D plot labels
    ttk.Label(label_window, text="1D Plot Title").grid(row=4, column=0)
    ttk.Entry(label_window, textvariable=label_1D_title).grid(row=4, column=1)

    ttk.Label(label_window, text="1D X Label").grid(row=5, column=0)
    ttk.Entry(label_window, textvariable=label_1D_xlabel).grid(row=5, column=1)

    ttk.Label(label_window, text="1D Y Label").grid(row=6, column=0)
    ttk.Entry(label_window, textvariable=label_1D_ylabel).grid(row=6, column=1)

    ttk.Label(label_window, text="1D Δx Series Label").grid(row=7, column=0)
    ttk.Entry(label_window, textvariable=label_1D_dx_label).grid(row=7, column=1)

    ttk.Label(label_window, text="1D Δy Series Label").grid(row=8, column=0)
    ttk.Entry(label_window, textvariable=label_1D_dy_label).grid(row=8, column=1)

    

def generate_IL_vs_dz():
    def plot_IL_vs_dz():
        try:
            Dx1 = float(entry_Dx1.get())
            Dx2 = float(entry_Dx2.get())
            Dy1 = float(entry_Dy1.get())
            Dy2 = float(entry_Dy2.get())
            wavelength = float(entry_wavelength.get())

            z_max = float(entry_zmax.get())
            z_step = float(entry_zstep.get())

            delta_z_values = np.arange(0, z_max + z_step, z_step)
            il_vs_z = []

            for dz in delta_z_values:
                eta_total = eta_delta_x(0, Dx1, Dx2) * eta_delta_y(0, Dy1, Dy2) * eta_ZM(Dx1, Dx2, Dy1, Dy2, dz, wavelength)
                il_vs_z.append(to_dB(eta_total))

            plt.figure()
            plt.plot(delta_z_values, il_vs_z, color='purple')
            plt.xlabel("Δz (μm)")
            plt.ylabel("Insertion Loss (dB)")
            #plt.title("Insertion Loss vs Axial Misalignment (Δz)")
            plt.axhline(-1, color='gray', linestyle='--', linewidth=0.8)
            plt.axhline(-3, color='gray', linestyle='--', linewidth=0.8)
            plt.grid(True)
            plt.show()
            z_window.destroy()

        except ValueError:
            print("Invalid input for Δz plot")

    # Create a small popup window
    z_window = tk.Toplevel(root)
    z_window.title("Δz Sweep Parameters")

    ttk.Label(z_window, text="Max Δz (μm):").grid(row=0, column=0)
    entry_zmax = ttk.Entry(z_window)
    entry_zmax.insert(0, "10")
    entry_zmax.grid(row=0, column=1)

    ttk.Label(z_window, text="Step (μm):").grid(row=1, column=0)
    entry_zstep = ttk.Entry(z_window)
    entry_zstep.insert(0, "0.1")
    entry_zstep.grid(row=1, column=1)

    ttk.Button(z_window, text="Generate Plot", command=plot_IL_vs_dz).grid(row=2, columnspan=2)



# Create the main window
root = tk.Tk()
root.title("Optical Coupling Calculator")

# Create and place the input fields and labels
ttk.Label(root, text="Dx1 (μm):").grid(row=0, column=0)
entry_Dx1_var = tk.StringVar(value="10")
entry_Dx1 = ttk.Entry(root, textvariable=entry_Dx1_var)
entry_Dx1.grid(row=0, column=1)

var_symmetry_Dx = tk.BooleanVar()
check_symmetry_Dx = ttk.Checkbutton(root, text="Symmetric", variable=var_symmetry_Dx, command=toggle_symmetry_Dx)
check_symmetry_Dx.grid(row=0, column=2)

ttk.Label(root, text="Dx2 (μm):").grid(row=1, column=0)
entry_Dx2_var = tk.StringVar(value="10")
entry_Dx2 = ttk.Entry(root, textvariable=entry_Dx2_var)
entry_Dx2.grid(row=1, column=1)

entry_Dx1_var.trace_add('write', sync_Dy1)
entry_Dx2_var.trace_add('write', sync_Dy2)


var_symmetry_Dy = tk.BooleanVar()
check_symmetry_Dy = ttk.Checkbutton(root, text="Symmetric", variable=var_symmetry_Dy, command=toggle_symmetry_Dy)
check_symmetry_Dy.grid(row=1, column=2)

ttk.Label(root, text="Dy1 (μm):").grid(row=2, column=0)
entry_Dy1_var = tk.StringVar(value="10")
entry_Dy1 = ttk.Entry(root, textvariable=entry_Dy1_var)
entry_Dy1.grid(row=2, column=1)

ttk.Label(root, text="Dy2 (μm):").grid(row=3, column=0)
entry_Dy2_var = tk.StringVar(value="10")
entry_Dy2 = ttk.Entry(root, textvariable=entry_Dy2_var)
entry_Dy2.grid(row=3, column=1)

ttk.Label(root, text="delta_x (μm):").grid(row=4, column=0)
entry_delta_x_var = tk.StringVar(value="5")
entry_delta_x = ttk.Entry(root, textvariable=entry_delta_x_var)
entry_delta_x.grid(row=4, column=1)

ttk.Label(root, text="delta_y (μm):").grid(row=5, column=0)
entry_delta_y_var = tk.StringVar(value="5")
entry_delta_y = ttk.Entry(root, textvariable=entry_delta_y_var)
entry_delta_y.grid(row=5, column=1)


ttk.Label(root, text="delta_teta (deg):").grid(row=6, column=0)
entry_delta_teta_var = tk.StringVar(value="0")
entry_delta_teta = ttk.Entry(root, textvariable=entry_delta_teta_var)
entry_delta_teta.grid(row=6, column=1)


ttk.Label(root, text="delta_phi (deg):").grid(row=7, column=0)
entry_delta_phi_var = tk.StringVar(value="0")
entry_delta_phi = ttk.Entry(root, textvariable=entry_delta_phi_var)
entry_delta_phi.grid(row=7, column=1)

ttk.Label(root, text="n (refractive index):").grid(row=8, column=0)
entry_n_var = tk.StringVar(value="1")
entry_n = ttk.Entry(root, textvariable=entry_n_var)
entry_n.grid(row=8, column=1)


ttk.Label(root, text="delta_z (μm):").grid(row=9, column=0)
entry_delta_z_var = tk.StringVar(value="0")
entry_delta_z = ttk.Entry(root, textvariable=entry_delta_z_var)
entry_delta_z.grid(row=9, column=1)


ttk.Label(root, text="Wavelength (μm):").grid(row=10, column=0)
entry_wavelength_var = tk.StringVar(value="1.55")
entry_wavelength = ttk.Entry(root, textvariable=entry_wavelength_var)
entry_wavelength.grid(row=10, column=1)

# Create and place the result labels
ttk.Label(root, text="η_M:").grid(row=11, column=0)
result_eta_M = tk.StringVar()
ttk.Label(root, textvariable=result_eta_M).grid(row=11, column=1)

ttk.Label(root, text="η_δx:").grid(row=12, column=0)
result_eta_delta_x = tk.StringVar()
ttk.Label(root, textvariable=result_eta_delta_x).grid(row=12, column=1)

ttk.Label(root, text="η_δy:").grid(row=13, column=0)
result_eta_delta_y = tk.StringVar()
ttk.Label(root, textvariable=result_eta_delta_y).grid(row=13, column=1)

ttk.Label(root, text="η_δteta:").grid(row=14, column=0)
result_eta_delta_teta = tk.StringVar()
ttk.Label(root, textvariable=result_eta_delta_teta).grid(row=14, column=1)

ttk.Label(root, text="η_δphi:").grid(row=15, column=0)
result_eta_delta_phi = tk.StringVar()
ttk.Label(root, textvariable=result_eta_delta_phi).grid(row=15, column=1)

ttk.Label(root, text="η_ZM:").grid(row=16, column=0)
result_eta_ZM = tk.StringVar()
ttk.Label(root, textvariable=result_eta_ZM).grid(row=16, column=1)

ttk.Label(root, text="η_total:").grid(row=17, column=0)
result_eta_total = tk.StringVar()
ttk.Label(root, textvariable=result_eta_total).grid(row=17, column=1)

# Create and place the calculate button
calculate_button = ttk.Button(root, text="Calculate", command=calculate_results)
calculate_button.grid(row=18, columnspan=2)

# Create and place the dB checkbox
var_dB = tk.BooleanVar()
check_dB = ttk.Checkbutton(root, text="Convert to dB", variable=var_dB, command=update_results)
check_dB.grid(row=19, columnspan=2)

# Create and place the heatmap button
ttk.Button(root, text="Generate IL map (2D)", command=generate_heatmap).grid(row=20, column=0)
ttk.Button(root, text="Generate IL map (1D)", command=generate_1D_IL_map).grid(row=20, column=1)

#label button
ttk.Button(root, text="Labels", command=open_label_settings).grid(row=20, column=2)

# IL vs Δz
ttk.Button(root, text="Generate IL vs Δz", command=generate_IL_vs_dz).grid(row=20, column=0, columnspan=3)



# Default label values
label_2D_title = tk.StringVar(value="Optical coupling efficiency")
label_2D_xlabel = tk.StringVar(value="δx (μm)")
label_2D_ylabel = tk.StringVar(value="δy (μm)")

label_1D_title = tk.StringVar(value="1D Insertion Loss vs Misalignment")
label_1D_xlabel = tk.StringVar(value="Δ (μm)")
label_1D_ylabel = tk.StringVar(value="Insertion Loss (dB)")
label_1D_dx_label = tk.StringVar(value="Δx Sweep")
label_1D_dy_label = tk.StringVar(value="Δy Sweep")



# Run the application
root.mainloop()
