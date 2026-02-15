# Intent



## Objective

- Reference: Asai Asaithambi, Numerical solution of the Burgers’ equation by automatic differentiation,
  Applied Mathematics and Computation, 216 (2010), 2700–2708.

- The objective of this study is to reinterpret the Taylor-series-based time integration method,
  originally proposed for time-dependent nonlinear PDEs, as a *pseudo-time relaxation and smoothing
  strategy* for elliptic problems.

- Specifically, we investigate the numerical properties of solving the Poisson equation by introducing a pseudo-time term and applying Taylor-series-based time integration.



## Goal:

- To clarify the fundamental convergence characteristics of the Taylor-series-based pseudo-time
  method for the Poisson equation, including its residual decay behavior and spectral smoothing
  properties.
- To evaluate the effectiveness of the Taylor method as a smoother within a multigrid framework,
  and to compare its performance with conventional smoothers such as SOR and SSOR.

- To propose and assess a *hierarchical Taylor strategy*, in which the Taylor expansion order and
  pseudo-time step size are adapted to the grid level (high-order/small Δt on fine grids,
  low-order/large Δt on coarse grids).

- To formulate and test a *Correction-Taylor* approach, where the multigrid correction equation
  is solved via pseudo-time Taylor integration, thereby unifying smoothing and coarse-grid correction within a common Taylor-based framework.





## Success Criteria:

- Compared with baseline methods, the proposed Taylor-based approaches achieve
  faster convergence in terms of iteration counts and/or wall-clock time.
- In multigrid configurations, the Taylor smoother demonstrates effective damping
  of high-frequency error components.

- The overall convergence behavior is robust and stable, with residuals exhibiting
  consistent long-term decay, even if strict monotonicity is not always observed.





## Comparison Baseline:

- Baseline solvers without multigrid:
  SOR, SSOR, and Conjugate Gradient (CG).
- Comparison metrics:
  - Residual history (L2 norm of residual)
  - Number of iterations to reach a prescribed tolerance
  - Wall-clock time
  - Error distribution compared to the exact solution

- Multigrid configurations:
  - Smoothers: SSOR vs. Taylor-based smoother
  - Fixed multigrid parameters vs. level-dependent Taylor parameters

- For CG, comparisons focus on overall computational efficiency rather than iteration counts,
  acknowledging the difference between direct elliptic solvers and pseudo-time relaxation methods.