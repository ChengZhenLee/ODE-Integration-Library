#pragma once
#include <vector>

namespace odelib {

template <typename Stepper, typename System, typename State, typename Time>
std::vector<State> integrate(Stepper stepper, System f, State y0, Time t0, Time t1, Time h) {
    std::vector<State> results{y0};
    State y = y0;
    Time t = t0;
    while (t < t1) {
        y = stepper.step(f, y, t, h);
        t += h;
        results.push_back(y);
    }
    return results;
};

}