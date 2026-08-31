#pragma once
#include <cmath>
#include <vector>

namespace odelib {

template <typename Stepper, typename System, typename State, typename Time>
std::vector<State> integrate(const Stepper& stepper, System&& f, State y0, Time t0, Time t1, Time h) {
    std::vector<State> res = {y0};

    long n = static_cast<long>(std::round((t1 - t0) / h));
    for (long i = 0; i < n; ++i) {
        Time t = t0 + static_cast<Time>(i) * h;
        res.push_back(stepper.step(f, res.back(), t, h));
    }

    return res;
};

}