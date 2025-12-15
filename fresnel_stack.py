import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import numpy as np

# --- 1. DATABASES ---

MATERIAL_DB = {
    "Air": 1.0,
    "Silica (SiO2)": 1.444,
    "Silicon (Si)": 3.47,
    "InP": 3.17,
    "GaAs": 3.3,
    "Germanium": 4.0,
    "BK7 Glass": 1.517,
    "Water": 1.33,
    "3D Polymer": 1.50,
    "SU-8": 1.59,
    "Si3N4": 2.0
}

# Consistent colors for materials
MATERIAL_COLORS = {
    "Air": "#ffffff",
    "Silica (SiO2)": "#4daf4a",      # Green
    "Silicon (Si)": "#0f5b78",       # Dark Blue
    "InP": "#f37735",                # Orange
    "GaAs": "#d62728",               # Red
    "Germanium": "#9467bd",          # Purple
    "BK7 Glass": "#aec7e8",          # Light Blue
    "Water": "#98df8a",              # Light Green
    "3D Polymer": "#17becf",         # Cyan
    "SU-8": "#e377c2",               # Pink
    "Si3N4": "#bcbd22"               # Olive
}
DEFAULT_COLOR = "#cccccc" # Gray for unknown materials

class FresnelStackTab(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        self.layers = []
        self.counter = 0 
        self.selected_layer_id = None # Track which layer is selected
        
        self.setup_ui()
        
        # Start with default
        self.add_layer("Air") 
        self.select_layer(self.layers[0]['id']) # Select the first one
        
        # --- BINDINGS ---
        # Bind keyboard events to the whole widget
        # Note: Widget needs focus to catch keys. We set focus when clicking layers.
        self.bind_all("<Up>", self.on_key_up)
        self.bind_all("<Down>", self.on_key_down)

    def setup_ui(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        
        # --- LEFT PANEL ---
        left_panel = ttk.Frame(self, padding=10)
        left_panel.grid(row=0, column=0, sticky="nsew")
        
        ttk.Label(left_panel, text="Layer Stack Builder", font=('Arial', 12, 'bold')).pack(pady=(0, 10), anchor='w')
        ttk.Label(left_panel, text="(Click graph to select, Arrow keys to move)", font=('Arial', 9, 'italic')).pack(pady=(0, 5), anchor='w')
        
        # Scrollable Area
        canvas_frame = ttk.Frame(left_panel)
        canvas_frame.pack(fill='both', expand=True, pady=5)
        
        self.canvas = tk.Canvas(canvas_frame, height=300, width=350)
        scrollbar = ttk.Scrollbar(canvas_frame, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Add Button
        btn_frame = ttk.Frame(left_panel)
        btn_frame.pack(fill='x', pady=5)
        ttk.Button(btn_frame, text="+ Add Layer", command=lambda: self.add_layer()).pack(side='left', fill='x', expand=True)

        # Results
        result_frame = ttk.LabelFrame(left_panel, text="Results", padding=10)
        result_frame.pack(fill='x', pady=10)
        
        self.lbl_trans = ttk.Label(result_frame, text="Transmission: -", font=('Arial', 10))
        self.lbl_trans.pack(anchor='w')
        
        self.lbl_loss = ttk.Label(result_frame, text="Total Loss: -", font=('Arial', 11, 'bold'))
        self.lbl_loss.pack(anchor='w', pady=(5,0))
        
        # --- RIGHT PANEL (Viz) ---
        right_panel = ttk.Frame(self, padding=10)
        right_panel.grid(row=0, column=1, sticky="nsew")
        
        self.fig, self.ax = plt.subplots(figsize=(4, 5))
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas_plot.get_tk_widget().pack(fill='both', expand=True)
        
        # Connect matplotlib click event
        self.canvas_plot.mpl_connect('button_press_event', self.on_plot_click)

    def get_color(self, material_name):
        return MATERIAL_COLORS.get(material_name, DEFAULT_COLOR)

    def add_layer(self, pre_selected_material=None):
        uid = self.counter
        self.counter += 1
        
        row_frame = ttk.Frame(self.scrollable_frame)
        row_frame.pack(fill='x', pady=2)
        
        # Make the row clickable to select
        row_frame.bind("<Button-1>", lambda e: self.select_layer(uid))
        
        # 1. Color Indicator (will update dynamically)
        lbl_color = tk.Label(row_frame, bg=DEFAULT_COLOR, width=2)
        lbl_color.pack(side='left', padx=(0,5))
        lbl_color.bind("<Button-1>", lambda e: self.select_layer(uid)) # Bind click
        
        # 2. Material Dropdown
        material_names = list(MATERIAL_DB.keys())
        combo_mat = ttk.Combobox(row_frame, values=material_names, width=15)
        combo_mat.pack(side='left', padx=2)
        
        # 3. Refractive Index Entry
        ttk.Label(row_frame, text="n:").pack(side='left', padx=2)
        entry_n = ttk.Entry(row_frame, width=7)
        entry_n.pack(side='left', padx=2)
        
        # --- Events ---
        def on_material_select(event):
            selected = combo_mat.get()
            if selected in MATERIAL_DB:
                entry_n.delete(0, tk.END)
                entry_n.insert(0, str(MATERIAL_DB[selected]))
                # Update color indicator immediately
                new_color = self.get_color(selected)
                lbl_color.config(bg=new_color)
                self.refresh_visualization()
                self.calculate_stack()

        combo_mat.bind("<<ComboboxSelected>>", on_material_select)
        # Bind click on widgets to selection too
        combo_mat.bind("<Button-1>", lambda e: self.select_layer(uid))
        entry_n.bind("<Button-1>", lambda e: self.select_layer(uid))
        entry_n.bind("<Return>", lambda e: self.calculate_stack())

        # Set default values
        if pre_selected_material and pre_selected_material in MATERIAL_DB:
            combo_mat.set(pre_selected_material)
            entry_n.insert(0, str(MATERIAL_DB[pre_selected_material]))
            lbl_color.config(bg=self.get_color(pre_selected_material))
        else:
            combo_mat.set("Select...")
            entry_n.insert(0, "1.0")

        # 4. Remove Button
        btn_del = ttk.Button(row_frame, text="X", width=3, command=lambda: self.remove_layer(uid))
        btn_del.pack(side='left', padx=5)
        
        # Store data
        self.layers.append({
            'id': uid,
            'widget': row_frame,
            'combo_name': combo_mat,
            'entry_n': entry_n,
            'lbl_color': lbl_color
        })
        
        self.select_layer(uid) # Auto select new layer
        self.refresh_visualization()
        self.calculate_stack()

    def select_layer(self, uid):
        self.selected_layer_id = uid
        self.focus_set() # Ensure keyboard events are caught
        self.refresh_visualization() # To draw the highlight border

    def on_plot_click(self, event):
        if event.inaxes != self.ax: return
        
        # Determine which bar was clicked.
        # Bars are at y = 0, 1, 2...
        # Since we plot reversed(layers), the visual top is the first in list.
        # Visualization Logic:
        # Index 0 (Top in list) -> Plotted at Y = N-1
        # Index N-1 (Bottom in list) -> Plotted at Y = 0
        
        N = len(self.layers)
        clicked_y = round(event.ydata)
        
        if 0 <= clicked_y < N:
            # Convert Y coordinate back to list index
            # y = N - 1 - index  => index = N - 1 - y
            list_index = (N - 1) - clicked_y
            
            if 0 <= list_index < N:
                uid = self.layers[list_index]['id']
                self.select_layer(uid)

    def on_key_up(self, event):
        self.move_selection(-1)

    def on_key_down(self, event):
        self.move_selection(1)

    def move_selection(self, direction):
        if self.selected_layer_id is None: return
        
        # Find current index
        idx = -1
        for i, layer in enumerate(self.layers):
            if layer['id'] == self.selected_layer_id:
                idx = i
                break
        
        if idx == -1: return
        
        new_idx = idx + direction
        
        # Check bounds and Swap
        if 0 <= new_idx < len(self.layers):
            self.layers[idx], self.layers[new_idx] = self.layers[new_idx], self.layers[idx]
            
            # Reorder widgets visually in the left panel
            for layer in self.layers:
                layer['widget'].pack_forget()
            for layer in self.layers:
                layer['widget'].pack(fill='x', pady=2)
                
            self.refresh_visualization()
            self.calculate_stack()

    def remove_layer(self, uid):
        for i, layer in enumerate(self.layers):
            if layer['id'] == uid:
                layer['widget'].destroy()
                self.layers.pop(i)
                break
        
        # If we removed selected layer, select the last one
        if self.selected_layer_id == uid:
            if self.layers:
                self.selected_layer_id = self.layers[-1]['id']
            else:
                self.selected_layer_id = None
                
        self.refresh_visualization()
        self.calculate_stack()

    def get_layer_values(self):
        values = []
        for layer in self.layers:
            try:
                name = layer['combo_name'].get()
                n_val = float(layer['entry_n'].get())
                values.append({'name': name, 'n': n_val, 'id': layer['id']})
            except ValueError:
                pass 
        return values

    def calculate_stack(self):
        layers = self.get_layer_values()
        if len(layers) < 2:
            self.lbl_trans.config(text="Total Transmission: -")
            self.lbl_loss.config(text="Total Loss: -")
            return

        total_transmission = 1.0
        
        for i in range(len(layers) - 1):
            n1 = layers[i]['n']
            n2 = layers[i+1]['n']
            
            if n1 + n2 == 0: R = 0
            else: R = ((n1 - n2) / (n1 + n2)) ** 2
            
            T = 1 - R
            total_transmission *= T
            
        if total_transmission <= 0:
            loss_db = 999
        else:
            loss_db = -10 * np.log10(total_transmission)
            
        self.lbl_trans.config(text=f"Total Transmission: {total_transmission:.4f}")
        self.lbl_loss.config(text=f"Total Loss: {loss_db:.4f} dB")

    def refresh_visualization(self):
        self.ax.clear()
        layers = self.get_layer_values()
        
        # Iterate reversed so Top Layer is at Top of Graph
        for i, layer in enumerate(reversed(layers)):
            name = layer['name']
            color = self.get_color(name)
            
            # Check if this layer is selected
            is_selected = (layer['id'] == self.selected_layer_id)
            edge_color = 'red' if is_selected else 'black'
            line_width = 3 if is_selected else 1
            
            self.ax.barh(y=i, width=1, height=1, color=color, edgecolor=edge_color, linewidth=line_width, align='center')
            
            # Text contrast logic
            text_color = 'white' if color in ['#0f5b78', '#d62728', '#9467bd'] else 'black'
            
            self.ax.text(0.5, i, f"{layer['name']}\nn={layer['n']}", 
                         ha='center', va='center', fontsize=9, fontweight='bold',
                         color=text_color)
            
        self.ax.set_ylim(-0.5, len(layers)-0.5)
        self.ax.set_xlim(0, 1)
        self.ax.axis('off')
        self.ax.set_title("Stack Cross-Section", fontsize=10)
        self.canvas_plot.draw()

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Fresnel Stack Test")
    root.geometry("800x500")
    tab = FresnelStackTab(root)
    tab.pack(fill='both', expand=True)
    root.mainloop()