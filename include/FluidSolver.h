#ifndef FLUIDSOLVER_H
#define FLUIDSOLVER_H

class FluidSolver
{
public:
    const static double DEFAULT_TIMESTEP;

    FluidSolver()
    {}
    FluidSolver(double timestep)
        : timestep(timestep)
    {}
    virtual ~FluidSolver()
    {}
    void update(int frameNumber);

protected:
    virtual void onAct(double time_step) = 0;
private:
    int frameIdx = -1;

    void act(double time_step);

    double timestep = DEFAULT_TIMESTEP;
};


#endif //FLUIDSOLVER_H
