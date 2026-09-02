#pragma once


namespace odelib {

template <typename State, typename Time>
struct RKDPStepper {
    template <typename System>
    State step(System&& f, const State& y, Time t, Time h) const {
        State k1 = f(t, y);
        State k2 = f(t + h/5.0, y + h/5.0 * k1);
        State k3 = f(t + h * 3.0/10.0, y + h * (3.0/40.0 * k1 + 9.0/40.0 * k2));
        State k4 = f(t + h * 4.0/5.0, y + h * (44.0/45.0 * k1 - 56.0/15.0 * k2 + 32.0/9.0 * k3));
        State k5 = f(t + h * 8.0/9.0, y + h * (19372.0/6561.0 * k1 - 25360.0/2187.0 * k2 + 64448.0/6561.0 * k3 - 212.0/729.0 * k4));
        State k6 = f(t + h, y + h * (9017.0/3168.0 * k1 - 355.0/33.0 * k2 + 46732.0/5247.0 * k3 + 49.0/176.0 * k4 - 5103.0/18656.0 * k5));

        State y5 = y + h * (35.0/384.0 * k1 + 500.0/1113.0 * k3 + 125.0/192.0 * k4 - 2187.0/6784.0 * k5 + 11.0/84.0 * k6);

        return y5;
    }
};

}