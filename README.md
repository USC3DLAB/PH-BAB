# phbab
## Progressive-Hedging Based Branch-and-Bound Algorithm (PH-BAB)

**Cite as:** Atakan, S., & Sen, S. (2018). A Progressive Hedging Based Branch-and-Bound Algorithm for Mixed-Integer Stochastic Programs. Available via Optimization Online.

**Introduction:** The PH-BAB algorithm can address multi-stage stochastic mixed-integer convex programs. The algorithm is described in Atakan & Sen (2018). This repository is the implementation of the algorithms described in Section 4 of that paper. The code is written in C++, and depends on CPLEX & Concert Technology. Further assumptions (some of which are not necessary, but speeded up the coding of the algorithm) are listed below.   

**Assumptions:**
- Solves two-stage stochastic mixed-integer programs
- Requires SMPS input
- 1st-stage variables must be pure binaries
- 1st-stage variables cannot have quadratic objective (can be easily relaxed)
- The convexity of the objective function and feasible sets is to the extent that CPLEX can handle
- Stochasticity appears only on the right-hand-side (can be easily relaxed)

**Parallelism:** This code can create multiple solvers for parallel-computation of subproblems. The code DOES NOT simultaneously solve multiple branch-and-bound nodes. Parallelization is achieved using Boost C++ libraries (Boost). For serial implementation (i.e., compiling w/o Boost Libraries) comment out "#define PHBB_parallel" in "config.h". In parallel mode, the code MAY NOT PRODUCE DETERMINISTIC RESULT. The reason is, subproblems are assigned to the first-available thread, and each optimization is warm-started with the basis from an earlier optimization (corresponding to an arbitrary subproblem). As a result, multiple optimal solutions / loose numerical tolerances / ... for nodal-relaxations can all lead to different outcomes in the branch-and-bound process.

If parallelism is desired, the following boost libraries must be linked to the project:
- libboost_atomic.a
- libboost_chrono.a
- libboost_system.a
- libboost_thread.a
- libboost_timer.a

**Notes:** The provided code has been tested on MacOS Sierra (v10.12.4) with CPLEX 12.6, and Boost v1.55.0.

**Example Usage:**

Solve sslp_5_50_100 using the PH-BAB algorithm (exact and default choice):
```
./ph-bb -f sslp_5_50_100
```

Solve sslp_5_50_100 using the PH-BAB-Apx algorithm (approximate alternative):
```
./ph-bb -f sslp_5_50_100 -a phbab-apx
```

Solve sslp_5_50_100, and simulate optimal 1st-stage solution using SSLP_5_50_2000.sto:
```
./ph-bb -f sslp_5_50_100 -s sslp_5_50_2000
```

Solve the entire SSLP dataset:
```
for f in ../datasets/sslp/*.cor; do
  probname=${f%.cor};
  ./ph-bb -f "$probname";
done
```
