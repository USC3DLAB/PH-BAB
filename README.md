# phbab
Progressive-Hedging Based Branch-and-Bound Algorithm (PH-BAB)

Cite as: Atakan, S., & Sen, S.. A Progressive Hedging Based Branch-and-Bound Algorithm for Mixed-Integer Stochastic Programs. Available via Optimization Online.

Introduction: The PH-BAB algorithm can address multi-stage stochastic mixed-integer convex programs. The algorithm is described in Atakan & Sen. This repository is the implementation of the algorithms described in Section 4 of that paper. Further assumptions (some of which are not necessary, but speeded up the coding of the algorithm) are listed below.   

Assumptions:
- Solves two-stage stochastic mixed-integer programs
- Requires SMPS input
- 1st-stage variables must be pure binaries
- 1st-stage variables cannot have quadratic objective (can be easily relaxed)
- The convexity of the objective function and feasible sets is to the extent that CPLEX 12.6 can handle
- Stochasticity appears only on the right-hand-side (can be easily relaxed)

Parallelism: This code can create multiple solvers for parallel-computation of subproblems. The code DOES NOT simultaneously solve multiple branch-and-bound nodes. Parallelization is achieved using Boost C++ libraries (Boost). For serial implementation (i.e., compiling w/o Boost Libraries) comment out "#define PHBB_parallel" in "config.h". In parallel mode, the code MAY NOT PRODUCE DETERMINISTIC RESULT. Subproblems are assigned to the first-available thread, therefore each optimization may be warm-started with the solution from an earlier computation. As a result, multiple optimal solutions / loose numerical tolerances for nodal-relaxations can all lead to different outcomes in the branch-and-bound process.

If parallelism is desired, the following boost libraries must be linked to the project:
- libboost_atomic.a
- libboost_chrono.a
- libboost_system.a
- libboost_thread.a
- libboost_timer.a

The provided code has been tested on MacOS Sierra (v10.12.4) with CPLEX 12.6, and Boost (v1.55.0).
