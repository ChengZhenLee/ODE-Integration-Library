# ODE Integration Library

A small, header-only C++20 library for numerically integrating ordinary differential equations (ODEs) and stochastic differential equations (SDEs).

## Features

- **Deterministic steppers**: Euler, RK4, Dormand-Prince (RKDP)
- **Adaptive step-size steppers**: Adaptive RK4, Adaptive RKDP
- **Stochastic steppers**: Euler-Maruyama, Milstein
- **Two-factor stochastic volatility**: Heston model (correlated Euler-Maruyama)
- Works with scalar states (`double`) or vector states (via [Eigen](https://eigen.tuxfamily.org/))
- Header-only — just include and go

## Layout

```
include/odelib/
  integrate.hpp        # generic integrate() driver functions
  steppers/             # individual stepper implementations
examples/               # small runnable examples for each stepper
tests/                  # unit tests
benchmarks/              # convergence/performance benchmarks + plotting scripts
```

## Requirements

- A C++20 compiler
- CMake 3.16+
- Python (only needed to regenerate benchmark plots)

## Building

```bash
cmake -B build
cmake --build build
```

This builds the examples, tests, and benchmarks.

## Running tests

```bash
ctest --test-dir build
```

## Running benchmarks

```bash
cmake --build build --target run_benchmarks
```

This runs all benchmark executables and regenerates the plots in `benchmarks/plots/`.

## Usage

Include a stepper and the `integrate` driver, then call `odelib::integrate`:

```cpp
#include "odelib/steppers/rk4.hpp"
#include "odelib/integrate.hpp"

auto decay = [](double t, double y) { return -y; };

int main() {
    double y0 = 1.0, t0 = 0.0, t1 = 5.0, h = 0.01;

    odelib::RK4Stepper<double> stepper;
    auto ys = odelib::integrate(stepper, decay, y0, t0, t1, h);
}
```

Adaptive steppers return a struct with `states`, `times`, and step sizes (`hs`) instead of a plain vector. Stochastic steppers (Euler-Maruyama, Milstein) take a drift and diffusion function in addition to the system function. The Heston stepper takes separate drift/diffusion functions for the asset price and variance processes, plus a correlation coefficient between the two Brownian motions. See the `examples/` directory for complete, runnable versions of each case (decay, harmonic oscillator, GBM, Black-Scholes, Heston).
