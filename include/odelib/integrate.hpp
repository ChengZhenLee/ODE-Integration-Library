#pragma once
#include <cmath>
#include <vector>


namespace odelib {

template <typename State, typename Time>
struct AdaptiveIntegrationResult {
    std::vector<State> states;
    std::vector<Time> times;
    std::vector<Time> hs;
};

template <typename Stepper, typename System, typename State, typename Time>
auto integrate(const Stepper& stepper, System&& f, State y0, Time t0, Time t1, Time h) {
    using StepReturnType = decltype(stepper.step(f, y0, t0, h));

    if constexpr(std::is_same_v<StepReturnType, State>) {
        std::vector<State> res = {y0};
        std::vector<Time> times = {t0};
        std::vector<Time> hs = {h};

        Time t = t0;

        while (t < t1) {
            res.push_back(stepper.step(f, res.back(), t, h));
            t += h;
        }

        return res;

    } else {
        std::vector<State> res = {y0};
        std::vector<Time> times = {t0};
        std::vector<Time> hs = {h};

        Time t = t0;
        Time hAdaptive = h;

        while (t < t1) {
            Time hClamped = std::min(hAdaptive, t1 - t);
            auto stepResult = stepper.step(f, res.back(), t, hClamped);

            if (stepResult.accepted) {
                t += stepResult.hUsed;
                res.push_back(stepResult.yNext);
                times.push_back(t);
                hs.push_back(stepResult.hUsed);
            }
            hAdaptive = stepResult.hSuggested;
        }

        return AdaptiveIntegrationResult<State, Time>{res, times, hs};
    }

};


template <typename Stepper, typename Drift, typename Diffusion, typename State, typename Time>
auto integrate(Stepper& stepper, Drift&& f, Diffusion&& g, State y0, Time t0, Time t1, Time h) {
    std::vector<State> res{y0};

    Time t = t0;

    while (t < t1) {
        res.push_back(stepper.step(f, g, res.back(), t, h));
        t += h;
    }

    return res;
};

}