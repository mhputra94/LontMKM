# cli.py
import sys
from .io_utils import parse_yaml, build_from_yaml, plot_solution, save_csv
from .core import solve_model

def main():
    if len(sys.argv) < 2:
        print("Usage: lontmkm <input.yaml>")
        sys.exit(1)

    data = parse_yaml(sys.argv[1])
    model, y0, (t0, tmax, npoints), prefix = build_from_yaml(data)
    sol = solve_model(model, y0, (t0, tmax), npoints)

    if sol.success:
        plot_solution(sol, model.species, f"{prefix}.png")
        save_csv(sol, model.species, f"{prefix}.csv")
        print(f"Results saved to {prefix}.png and {prefix}.csv")
    else:
        print("Integration failed:", sol.message)

if __name__ == "__main__":
    main()
