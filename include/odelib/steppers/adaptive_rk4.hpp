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
            Time hSuggested;
            bool accepted;
            
            // Expand the step-size if the error magnitude is a wide margin lower than the tolerance
            if (errorMagnitude < tol / 10.0) {
                hSuggested = h * 1.5;
                accepted = true;
            // Keep the suggested step size constant if not far below the tolerance
            } else if (errorMagnitude < tol) {
                hSuggested = h;
                accepted = true;
            // Suggest a smaller step size if the error exceeds tolerance
            } else {
                hSuggested = h / 2.0;
                accepted = false;
            }

            return { y_next, h, hSuggested, accepted };
        }
};

}