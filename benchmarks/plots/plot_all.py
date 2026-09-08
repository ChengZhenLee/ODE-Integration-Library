import os

import matplotlib.pyplot as plt
import pandas as pd

PLOTS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(PLOTS_DIR, "..", "data"))


def plot_lines(df, x_col, y_col, xlabel, ylabel, title, filename,
                xscale="linear", yscale="linear"):
    """Plot one line per `method` group from df, x_col vs y_col.

    Sorts each group by x_col before plotting, since CSV row order
    (e.g. sorted by tolerance) does not always match x-axis order
    (e.g. time_seconds, which is noisy) -- unsorted data produces a
    jagged, back-and-forth line instead of a clean curve.
    """
    fig, ax = plt.subplots()
    for method_name, group in df.groupby("method"):
        group = group.sort_values(x_col)
        ax.plot(group[x_col], group[y_col], marker='o', label=method_name)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    fig.savefig(os.path.join(PLOTS_DIR, filename))
    plt.close(fig)


def plot_fixed_step_convergence():
    df = pd.read_csv(f"{DATA_DIR}/FixedStepConvergence.csv")

    plot_lines(df, "h", "error",
               "Step size (h)", "Error", "Fixed-Step Convergence",
               "fixed_step_convergence.png",
               xscale="log", yscale="log")

    plot_lines(df, "time_seconds", "error",
               "Time (seconds)", "Error", "Work-Precision Diagram",
               "fixed_step_convergence_work_precision.png",
               xscale="log", yscale="log")


def plot_adaptive_step_trace():
    df = pd.read_csv(f"{DATA_DIR}/AdaptiveStepTrace.csv")

    plot_lines(df, "t", "h",
               "Time (t)", "Step size (h)", "Adaptive Step Trace",
               "adaptive_step_trace.png")


def plot_adaptive_step_tolerance():
    df = pd.read_csv(f"{DATA_DIR}/AdaptiveStepTolerance.csv")

    plot_lines(df, "tolerance", "total_steps",
               "Tolerance", "Total Steps", "Adaptive Step Tolerance",
               "adaptive_step_tolerance.png",
               xscale="log")

    plot_lines(df, "time_seconds", "error",
               "Time (seconds)", "Error", "Work-Precision Diagram (Adaptive)",
               "adaptive_step_tolerance_work_precision.png",
               xscale="log", yscale="log")


def plot_strong_convergence():
    df = pd.read_csv(f"{DATA_DIR}/StrongConvergencePlot.csv")

    plot_lines(df, "h", "error",
               "Step size (h)", "Strong Error",
               "EM vs Milstein Strong Convergence",
               "strong_convergence.png",
               xscale="log", yscale="log")


def plot_monte_carlo_n():
    df = pd.read_csv(f"{DATA_DIR}/MonteCarloN.csv")

    plot_lines(df, "N", "error",
               "Number of paths (N)", "Pricing Error",
               "Monte Carlo Convergence",
               "monte_carlo_convergence.png",
               xscale="log", yscale="log")


def plot_thread_speedup():
    df = pd.read_csv(f"{DATA_DIR}/ThreadSpeedup.csv")

    # Speedup = time at 1 thread / time at N threads, computed here rather
    # than stored in the CSV, since it's a derived quantity.
    df = df.sort_values("thread_count")
    baseline = df[df["thread_count"] == 1]["time_seconds"].iloc[0]
    df["speedup"] = baseline / df["time_seconds"]

    fig, ax = plt.subplots()
    ax.plot(df["thread_count"], df["speedup"], marker='o', label="Measured speedup")
    ax.plot(df["thread_count"], df["thread_count"], linestyle='--', label="Ideal (linear) speedup")

    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Speedup")
    ax.set_title("Multithreading Speedup")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.savefig(os.path.join(PLOTS_DIR, "thread_speedup.png"))
    plt.close(fig)


def plot_walk_comparison():
    df = pd.read_csv(f"{DATA_DIR}/WalkComparison.csv")

    plot_lines(df, "t", "state",
               "Time (t)", "Stock price (S)",
               "Euler-Maruyama vs Milstein: Single Path",
               "walk_comparison.png")


if __name__ == "__main__":
    plot_fixed_step_convergence()
    plot_adaptive_step_trace()
    plot_adaptive_step_tolerance()
    plot_strong_convergence()
    plot_monte_carlo_n()
    plot_thread_speedup()
    plot_walk_comparison()
    print("All plots generated.")