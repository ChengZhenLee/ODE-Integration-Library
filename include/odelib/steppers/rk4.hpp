#pragma once


namespace odelib {

template <typename State, typename Time=double>
struct RK4Stepper {
    template <typename System>
    State step(System&& f, const State& y, Time t, Time h) const {
        State k1 = f(t, y);
        State k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
        State k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
        State k4 = f(t + h, y + h * k3);

        return y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
};

}