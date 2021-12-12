I created this algorithm to optimize a Model Predictive Control (MPC) parameters, applied to Doubly Fed Induction Generator (DFIG).

Here some characteristics of my algorithm:

- The MPC and DFIG model are implemented in a Simulink Simulation (`dfig.slx`);

- The MPC parameters to be optimized are the weighing matrices;

- It uses a variant of Particle Swarm Optimization (PSO) algorithm, called wDOCHM-PSO (It supports constraints and adaptive inertia);
- You can **adapt my algorithm to any Simulink Simulation** , to do that, just modify `dfig.slx` and `fitness.m`.
- You can **change the optimization algorithm** to any other meta-heuristic or evolutionary algorithm, like Genetic Algorithm, by replacing or modifying `pso.m`.

For more details, or if you want to use my algorithm don't forget to cite my paper (It is open access):

> **L. L. Rodrigues, J. S. Solís-Chaves, O. A. C. Vilcanqui and A. J. S. Filho, "Predictive Incremental Vector Control for DFIG With Weighted-Dynamic Objective Constraint-Handling Method-PSO Weighting Matrices Design," in *IEEE Access*, vol. 8, pp. 114112-114122, 2020, doi: 10.1109/ACCESS.2020.3003285.**
> Abstract: This paper proposes a Particle Swarm Optimization (PSO) based method, the Weighted-Dynamic-Objective Constraint-Handling PSO Method (WDOCHM-PSO). This was used to design the weighting matrices of an incremental Model-Based Predictive Controller (MBPC) for a Doubly Fed Induction Generator (DFIG) applied in a small-scale wind energy system. In contrast to the original PSO, the proposed method has an inner mechanism for dealing with constraints and an adaptive search factor. Additionally, the proposed incremental MPBC implementation does not need the flux information, since the intrinsic integral action rejects the constant flux disturbance. Finally, experimental results show that the proposed controller with the new constraint handling design method is nearly two times faster (In terms of settling time) than other formulations reported in the literature.
> URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9121987&isnumber=8948470

**For more details** see my paper: https://ieeexplore.ieee.org/document/9121987


**TL; DR:** This algorithm uses PSO to optimize a Simulink Simulation parameters. If you want, you can use it to any simulink simulation, you just need to modify or replace `dfig.slx`

# Instructions
## Testing

If you want to test, follow the instructions bellow:


1. Save all files in the same folder;
1. Open `dfig.slx`;
2. Open `fitness.m`;
2. Run `fitness.m` and see the output;
2. Run `pso.m` to find optimal parameters (It would take several hours);

## Fitness parameters

`fitness.m`  parameters for MPC and DFIG simulation details are:

```matlab
Ny=20; %Predict horizon
Nu=10; %Control horizon
Variac=120; %DC Link voltage
StepInit=[1,1]; %Beginning current
StepFinal=[3,3]; %Final current
StepTime=[1,1.5]; %Step time (Direct and quadrature current)
```

`fitness.m` outputs are:

```matlab
f %fitness value
phi %Constraint fitness value
tss %Settling time
ITAE %Performance index
Mp %Overshoot
sucess %Simulation was successful
Signal %Currents inputs (For plotting)
```

## PSO parameters

`pso.m` parameters are:

```matlab
K = 16;   % Population
N  = 300; % Generations
cg = 1.5; % Global weighing
cp = 1.5; % Personal weighting
```
