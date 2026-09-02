#pragma once

#include <random>


namespace odelib {

template <typename State, typename Time=double>
struct EulerMaruyamaStepper { 
    private:
        std::mt19937_64 engine;

    public:
        EulerMaruyamaStepper() {}

        EulerMaruyamaStepper(unsigned int seed) : engine(seed) {}

        template <typename Drift, typename Diffusion>
        State step(Drift&& f, Diffusion&& g, State y, Time t, Time h) {
            // Random draw follows normal distribution with mean 0 and standard deviation sqrt(h)
            std::normal_distribution<Time> distribution(0, std::sqrt(h));

            return y + f(t, y) * h + g(t, y) * distribution(engine);
        }

        template <typename Drift, typename Diffusion>
        State step(Drift&& f, Diffusion&& g, State y, Time t, Time h, Time dW) {
            return y + f(t, y) * h + g(t, y) * dW;
        }
};

}