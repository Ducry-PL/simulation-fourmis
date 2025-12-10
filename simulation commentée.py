import random
from math import pow
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from matplotlib.widgets import CheckButtons
import tkinter as tk
from tkinter import messagebox

# ==========================================
# 1. FENÊTRE DE CONFIGURATION
# ==========================================

def lancer_configuration():
    """Ouvre une fenêtre pour paramétrer la simulation."""
    params = {}
    
    def valider():
        try:
            params['NCmax'] = int(entry_cycles.get())       # Nombre de cycles total
            params['alpha'] = float(entry_alpha.get())      # Exposant alpha (poids de la trace tau)
            params['beta'] = float(entry_beta.get())        # Exposant beta (poids de la visibilité eta)
            params['evaporation'] = float(entry_evap.get()) # Variable rho (taux d'évaporation)
            params['m'] = int(entry_ants.get())             # Variable m (nombre de fourmis)
            root.destroy()
        except ValueError:
            messagebox.showerror("Erreur", "Veuillez entrer des nombres valides.")

    root = tk.Tk()
    root.title("Configuration ACO")
    root.geometry("300x400")
    root.configure(bg='#222')

    lbl_style = {'bg': '#222', 'fg': 'white', 'font': ('Arial', 10)}
    title_style = {'bg': '#222', 'fg': 'white', 'font': ('Arial', 12, 'bold')}
    entry_bg = '#333'
    entry_fg = 'white'

    tk.Label(root, text="Paramètres de la Simulation", **title_style).pack(pady=15)

    # --- Champs ---
    tk.Label(root, text="Nombre de cycles (NCmax) :", **lbl_style).pack()
    entry_cycles = tk.Entry(root, bg=entry_bg, fg=entry_fg, insertbackground='white')
    entry_cycles.insert(0, "100")
    entry_cycles.pack(pady=5)

    tk.Label(root, text="Influence Phéromones (Alpha) :", **lbl_style).pack()
    entry_alpha = tk.Entry(root, bg=entry_bg, fg=entry_fg, insertbackground='white')
    entry_alpha.insert(0, "1.0")
    entry_alpha.pack(pady=5)

    tk.Label(root, text="Influence Visibilité (Beta) :", **lbl_style).pack()
    entry_beta = tk.Entry(root, bg=entry_bg, fg=entry_fg, insertbackground='white')
    entry_beta.insert(0, "5.0")
    entry_beta.pack(pady=5)

    tk.Label(root, text="Taux évaporation (0.0 à 1.0) :", **lbl_style).pack()
    entry_evap = tk.Entry(root, bg=entry_bg, fg=entry_fg, insertbackground='white')
    entry_evap.insert(0, "0.5")
    entry_evap.pack(pady=5)
    
    tk.Label(root, text="Nombre de fourmis :", **lbl_style).pack()
    entry_ants = tk.Entry(root, bg=entry_bg, fg=entry_fg, insertbackground='white')
    entry_ants.insert(0, "6")
    entry_ants.pack(pady=5)

    btn = tk.Button(root, text="LANCER LA SIMULATION", command=valider, 
                    bg='#00d2ff', fg='black', font=('Arial', 10, 'bold'))
    btn.pack(pady=20)

    root.mainloop()
    
    if not params: exit() 
    return params

user_params = lancer_configuration()

# ==========================================
# 2. INITIALISATION
# ==========================================

plt.style.use('dark_background')

selected_cities = ["Paris", "Brest", "Lyon", "Nice", "Pau", "Troyes"]
full_distances = {
  "Paris":  {"Paris":0,"Brest":7,"Lyon":7,"Nice":8,"Pau":13,"Troyes":1},
  "Brest":  {"Paris":7,"Brest":0,"Lyon":11,"Nice":15,"Pau":13,"Troyes":8},
  "Lyon":   {"Paris":7,"Brest":11,"Lyon":0,"Nice":4,"Pau":7,"Troyes":6},
  "Nice":   {"Paris":8,"Brest":15,"Lyon":4,"Nice":0,"Pau":7,"Troyes":9},
  "Pau":    {"Paris":13,"Brest":13,"Lyon":7,"Nice":7,"Pau":0,"Troyes":10},
  "Troyes": {"Paris":1,"Brest":8,"Lyon":6,"Nice":9,"Pau":10,"Troyes":0}
}
distances = {city: {target: full_distances[city][target] for target in selected_cities} for city in selected_cities}

city_coords = {
    "Brest":  (0, 10), "Lyon": (8.66, 5), "Nice": (8.66, -5),
    "Pau":    (0, -10), "Troyes": (-8.66, -5), "Paris":  (-8.66, 5)
}

cities = list(distances.keys())
n = len(cities) # n correspond au nombre de villes (sommets du graphe)

# Récupération des variables du problème
NCmax = user_params['NCmax']
alpha = user_params['alpha']        # Exposant pour l'importance de la phéromone
beta = user_params['beta']          # Exposant pour l'importance de la visibilité
evaporation = user_params['evaporation'] # Correspond à rho (facteur de diminution)
m = user_params['m']                # Nombre de fourmis par cycle

Q = 100.0   # Constante Q pour le dépôt (Delta Tau = Q/L)
tau0 = 0.1  # Valeur initiale des phéromones (Tau_0)

# Matrice des distances (d_ij)
dist = [[float(distances[ci][cj]) for cj in cities] for ci in cities]

# Calcul de la matrice de visibilité (Eta_ij = 1 / d_ij)
# On parcourt chaque paire de villes i, j pour remplir la matrice
eta = [[0.0 if i == j or dist[i][j] == 0 else 1.0 / dist[i][j] for j in range(n)] for i in range(n)]

# Matrice des phéromones (Tau_ij) initialisée à tau0
tau = [[tau0]*n for _ in range(n)]

global_best_len = float('inf')
global_best_tour = []
show_pheromones = True

# ==========================================
# 3. LOGIQUE ACO
# ==========================================

def roulette_choose(scores, candidates):
    """Sélectionne une ville selon sa probabilité (Roulette Wheel)."""
    # total correspond au dénominateur de la formule de probabilité (somme des scores)
    total = sum(scores[j] for j in candidates)
    if total == 0.0: return random.choice(candidates)
    r = random.random() * total
    cum = 0.0
    # On parcourt les candidats pour trouver où tombe r dans la somme cumulée
    for j in candidates:
        cum += scores[j]
        if r <= cum: return j
    return candidates[-1]

def build_tour(start_idx):
    """Construit un tour complet pour une fourmi."""
    visited = [False]*n
    tour = [start_idx]
    visited[start_idx] = True
    current = start_idx
    
    # Tant que la fourmi n'a pas visité les n villes
    while len(tour) < n:
        # Liste des villes non encore visitées (N_i^k)
        candidates = [j for j in range(n) if not visited[j]]
        scores = [0.0] * n
        # Boucle sur chaque ville candidate j pour calculer son attractivité
        for j in candidates:
            # Calcul du numérateur de la probabilité : [tau_ij]^alpha * [eta_ij]^beta
            scores[j] = pow(tau[current][j], alpha) * pow(eta[current][j], beta)
        
        next_j = roulette_choose(scores, candidates)
        tour.append(next_j)
        visited[next_j] = True
        current = next_j
    
    tour.append(start_idx) # La fourmi retourne au point de départ
    return tour

def tour_length(tour):
    """Calcule la longueur totale L_k d'un tour."""
    # Somme des distances d_ij pour chaque arête du tour
    return sum(dist[tour[k]][tour[k+1]] for k in range(len(tour)-1))

# ==========================================
# 4. VISUALISATION (Matplotlib)
# ==========================================

fig, ax = plt.subplots(figsize=(10, 10))
fig.canvas.manager.set_window_title(f"Simulation ACO")
fig.patch.set_facecolor('#121212')
ax.set_facecolor('#121212')

xs = [city_coords[c][0] for c in cities]
ys = [city_coords[c][1] for c in cities]
ax.set_xlim(-15, 15); ax.set_ylim(-15, 15)
ax.axis('off')

# Affichage des villes
ax.scatter(xs, ys, s=500, c='#00d2ff', zorder=10, edgecolors='white')
for i, city in enumerate(cities):
    ax.text(xs[i], ys[i]+1.8, city, color='white', ha='center', fontsize=11, fontweight='bold')

# Initialisation des lignes pour les traces de phéromones
edge_lines = {}
for i in range(n):
    for j in range(i+1, n):
        line, = ax.plot([], [], color='#2ecc71', alpha=0.0, linewidth=0.5)
        edge_lines[(i,j)] = line

best_tour_line, = ax.plot([], [], color='#ff0055', linewidth=4.0, zorder=5)

info_text = ax.text(0.02, 0.95, "", transform=ax.transAxes, color='white', 
                    bbox=dict(facecolor='#222', alpha=0.8))

# Checkbox
ax_check = plt.axes([0.05, 0.05, 0.25, 0.05], frameon=False)
check = CheckButtons(ax_check, ['Afficher Phéromones'], [True])
def toggle_pheromones(label):
    global show_pheromones
    show_pheromones = not show_pheromones
check.on_clicked(toggle_pheromones)

def update(frame):
    global global_best_len, global_best_tour, tau
    
    # Delta correspond à la matrice des dépôts de ce cycle (Delta Tau_ij)
    delta = [[0.0]*n for _ in range(n)]
    
    # Boucle sur chaque fourmi k (de 0 à m-1)
    for k in range(m):
        tour = build_tour(k % n)
        L = tour_length(tour) # L_k : longueur du tour de la fourmi k
        
        if L < global_best_len:
            global_best_len = L
            global_best_tour = tour[:]
            
        # Calcul de la quantité de phéromone à déposer : Delta Tau = Q / L_k
        deposit = Q / L
        # On parcourt les arêtes (u, v) empruntées par la fourmi k
        for e in range(len(tour)-1):
            u, v = tour[e], tour[e+1]
            delta[u][v] += deposit
            delta[v][u] += deposit

    # Mise à jour globale des phéromones sur toutes les arêtes
    max_tau_curr = 0
    for i in range(n):
        for j in range(n):
            # Application de l'évaporation (1-rho) et ajout des dépôts (delta)
            tau[i][j] = (1.0 - evaporation) * tau[i][j] + delta[i][j]
            max_tau_curr = max(max_tau_curr, tau[i][j])

    # Mise à jour graphique
    if max_tau_curr == 0: max_tau_curr = 1
    
    # On parcourt chaque ligne représentant une arête pour ajuster son épaisseur
    for (i, j), line in edge_lines.items():
        if not show_pheromones:
            line.set_alpha(0)
        else:
            ratio = tau[i][j] / max_tau_curr
            if ratio < 0.05:
                line.set_alpha(0)
            else:
                line.set_data([xs[i], xs[j]], [ys[i], ys[j]])
                line.set_linewidth(1.0 + 6 * ratio)
                line.set_alpha(0.1 + 0.9 * ratio)

    if global_best_tour:
        bx = [xs[idx] for idx in global_best_tour]
        by = [ys[idx] for idx in global_best_tour]
        best_tour_line.set_data(bx, by)

    info_text.set_text(f"CYCLE: {frame+1}/{NCmax}\nRecord: {int(global_best_len)}")
    return list(edge_lines.values()) + [best_tour_line, info_text]

ani = animation.FuncAnimation(fig, update, frames=NCmax, interval=100, blit=False, repeat=False)
plt.show()