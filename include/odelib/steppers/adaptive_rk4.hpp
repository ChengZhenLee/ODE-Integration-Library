#pragma once


#include "rk4.hpp"

namespace odelib {

template <typename State, typename Time>
struct AdaptiveRK4Stepper {
    private:
        double tol;

        double error_magnitude(const State& diff) const {
            if constexpr (requires { diff.norm(); }) {
                return diff.norm();
            } else {
                return std::abs(diff);
            }
        }

    public:
        AdaptiveRK4Stepper(double t = 1e-8) : tol(t) {};

        struct StepResult {
            State yNext;
            Time hUsed;
            Time hSuggested;
            bool accepted;
        };

        template <typename System>
        StepResult step(System&& f, const State& y, Time t, Time h) const {
            RK4Stepper<State, Time> stepper;
            State y_next = stepper.step(f, y, t, h);
            State y_smaller = stepper.step(f, stepper.step(f, y, t, h/2), t + h/2, h/2);

            double errorMagnitude = error_magnitude((y_next - y_smaller));
            Time h_next;
            bool accepted;
            
            if (errorMagnitude < tol / 10) {
                h_next = h * 1.5;
                accepted = true;
            } else if (errorMagnitude < tol) {
                h_next = h;
                accepted = true;
            } else {
                h_next = h/2;
                accepted = false;
            }

            return { y_next, h, h_next, accepted };
        }
};

}