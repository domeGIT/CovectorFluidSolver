#include <iostream>
#include "FluidSolver.h"

const double FluidSolver::DEFAULT_TIMESTEP = 0.01;
extern double subDt;

void FluidSolver::update(int frameNumber)
{
    
    int numFrames = 0;
    if (frameNumber > frameIdx)
    {
        numFrames = frameNumber - frameIdx;
        for (int i = frameIdx; i < frameIdx + numFrames; i++)
        {
            act(timestep);
            //timestep = subDt;
        }
        frameIdx = frameNumber;
    }

}

void FluidSolver::act(double time_step)
{
        onAct(time_step);
}
