#pragma once


namespace odelib {

template <typename State, typename Time=double>
struct EulerStepper {
    template <typename System>
    State step(System&& f, const State& y, Time t, Time h) const {
        return y + h * f(t, y);
    }
};

}