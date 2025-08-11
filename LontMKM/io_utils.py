# io_utils.py
import yaml
import matplotlib.pyplot as plt
import csv
from .core import MicrokineticModel, Reaction

def parse_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def build_from_yaml(data):
    species = data['species']
    global_T = float(data.get('global', {}).get('T', 298.0))
    reactions = [Reaction(r['name'], r['stoichiometry'], r['type'], r.get('params', {}), global_T) for r in data['reactions']]
    model = MicrokineticModel(species, reactions)
    y0 = [float(data.get('initial_conditions', {}).get(s, 0.0)) for s in species]
    t0 = float(data.get('time', {}).get('t0', 0.0))
    tmax = float(data.get('time', {}).get('tmax', 1.0))
    npoints = int(data.get('time', {}).get('npoints', 1000))
    output_prefix = data.get('output_prefix', 'microkinetic')
    return model, y0, (t0, tmax, npoints), output_prefix

def plot_solution(sol, species, filename):
    plt.figure(figsize=(8,5))
    for i, s in enumerate(species):
        plt.plot(sol.t, sol.y[i], label=s)
    plt.xlabel('time (s)')
    plt.ylabel('concentration')
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

def save_csv(sol, species, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["time"] + species)
        for i in range(len(sol.t)):
            writer.writerow([sol.t[i]] + [sol.y[j][i] for j in range(len(species))])

