template <typename State, template <typename> typename Stepper>
class ODEIntegrator : private Stepper<State>
{
    public:
    auto solve(State state, State (*derivative)(const State&), double start, double end) -> State
    {
        for (auto current = start; current < end; current += step_size)
        {
            if ((current + step_size - end)*(end - start) > 0.0) // if step overshoots
            {
                step_size = end - current;
            }
            state = this->step(state, derivative, step_size);
        }
        return state;
    }
    private:
    double step_size{1.0};
};

template <typename State>
struct EulerStepper
{
    auto step(State state, State (*derivative)(const State&), double& step_size) -> State
    {
        return state + derivative(state)*step_size;
    }
};

template <typename State>
struct RK2Stepper
{
    auto step(State state, State (*derivative)(const State&), double& step_size) -> State
    {
        auto k1 = derivative(state)*step_size;
        return state + derivative(state + 0.5*k1)*step_size;
    }
};

template <typename State>
struct RK4Stepper
{
    auto step(State state, State (*derivative)(const State&), double& step_size) -> State
    {
        auto k1 = derivative(state)*step_size;
        auto k2 = derivative(state + 0.5*k1)*step_size;
        auto k3 = derivative(state + 0.5*k2)*step_size;
        auto k4 = derivative(state + k3)*step_size;
        return state + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
    }
};