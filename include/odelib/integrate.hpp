#pragma once
#include <vector>

namespace odelib {

template <typename Stepper, typename System, typename State, typename Time>
std::vector<State> integrate(const Stepper& stepper, System&& f, State y0, Time t0, Time t1, Time h) {
    Time t = t0;
    std::vector<State> res = {y0};

    while (t <= t1) {
        res.push_back(stepper.step(f, res.back(), t, h));
        t += h;
    }

    return res;
};

}