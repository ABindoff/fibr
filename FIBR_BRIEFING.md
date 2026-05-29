# fibr: Holonomy Diagnostics for Hierarchical MCMC
## Project Briefing for Claude Code

---

## What This Project Is

`fibr` is an R package (with compiled backend to follow) that diagnoses a specific geometric pathology in MCMC chains for hierarchical models: **non-trivial holonomy** arising from the fiber bundle structure of the parameter space.

This is not a standard convergence diagnostic. It is detecting something that R-hat, ESS, and divergence counts cannot see — systematic rotation of local parameters as global parameters traverse closed loops in hyperparameter space.

---

## The Geometric Setup

A hierarchical model has a natural fiber bundle structure:

- **Base space** $B$: hyperparameters (e.g. $\mu, \sigma$ in a normal-normal hierarchy)
- **Fiber** $F_\theta$: local/group-level parameters conditioned on hyperparameters (e.g. $\alpha_1, \ldots, \alpha_J$)
- **Total space** $E = B \times F$: the full joint parameter space

The Fisher information metric on the total space has off-diagonal blocks coupling base and fiber directions. These off-diagonal blocks define a **connection** on the bundle — a notion of "horizontal" (base-driven) vs "vertical" (fiber-driven) motion.

The **curvature** of this connection means that parallel-transporting a fiber point around a closed loop in $B$ does not return it to the same position. This rotation is the **holonomy**, and it is a geometric obstruction to efficient MCMC mixing that standard samplers (including NUTS) do not account for.

The centering/non-centering reparameterisation debate in hierarchical models is, geometrically, a search for a **trivialisation** of the bundle that flattens the connection. It works when the model is well-identified; when it is not, no trivialisation is flat and you need something more principled.

---

## The Toy Model

A two-level logistic regression (GLMM):

$$y_{ij} \sim \text{Bernoulli}(\text{logit}^{-1}(\alpha_i + x_{ij}^\top \beta))$$
$$\alpha_i \sim \mathcal{N}(\mu, \sigma^2), \quad i = 1, \ldots, J$$

Parameter roles:
- **Base space**: $(\mu, \sigma)$ — the hyperparameters
- **Fiber**: $(\alpha_1, \ldots, \alpha_J)$ — group-level intercepts
- **Fixed effects**: $\beta$ (treat as part of the base for now, or marginalise)

Start with:
- $J = 8$ groups
- $n_j \approx 20$ observations per group
- Simulated data (see simulation spec below)

This model is ideal because:
1. The fiber geometry is non-trivial — $\alpha_i$ are correlated through the likelihood even though the prior factorises
2. The connection has an analytic form we can check the diagnostic against
3. It is a canonical use case where non-centering is known to matter

---

## First Milestone: Holonomy Diagnostic

### What It Measures

Given a posterior chain on $(\mu, \sigma, \alpha_1, \ldots, \alpha_J)$:

1. **Identify approximate loops in base space**: Find segments of the chain where $(\mu, \sigma)$ starts and ends near the same point (within some tolerance $\varepsilon$)
2. **Extract fiber displacement**: For each loop, record $\alpha_{\text{start}}$ and $\alpha_{\text{end}}$
3. **Estimate the transport map**: Fit $\alpha_{\text{end}} \approx H \cdot \alpha_{\text{start}}$ across many loops via regression or SVD
4. **Measure holonomy**: Decompose $H$ — if $H \approx I$, the connection is flat; eigenvalues away from 1 (especially complex pairs indicating rotation) are the holonomy signal

The output should be:
- A holonomy matrix $\hat{H}$ with uncertainty
- Its eigenspectrum (the key diagnostic)
- A scalar summary: e.g. $\|H - I\|_F$ or the angle of the dominant rotation

### Key Design Decisions

- Loop detection: use $k$-nearest neighbours in $(\mu, \sigma)$ space, not exact returns (chains never return exactly). Tolerance $\varepsilon$ is a tuning parameter.
- The transport map regression should be weighted by loop length / quality
- Need to separate genuine holonomy from noise — bootstrap the eigenvalues

---

## Simulation Spec

```r
set.seed(42)
J <- 8        # groups
n_j <- 20     # obs per group
mu_true <- 0
sigma_true <- 1.5
beta_true <- c(0.8, -0.5)  # two fixed predictors

# Simulate
alpha_true <- rnorm(J, mu_true, sigma_true)
group_id <- rep(1:J, each = n_j)
X <- matrix(rnorm(J * n_j * 2), ncol = 2)
eta <- alpha_true[group_id] + X %*% beta_true
y <- rbinom(J * n_j, 1, plogis(eta))
```

---

## Stan Model

Write a centred parameterisation first (non-centred comes later for comparison):

```stan
data {
  int<lower=0> N;
  int<lower=1> J;
  array[N] int<lower=1,upper=J> group;
  matrix[N, 2] X;
  array[N] int<lower=0,upper=1> y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector[J] alpha;
  vector[2] beta;
}
model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  alpha ~ normal(mu, sigma);
  beta ~ normal(0, 2);
  y ~ bernoulli_logit(alpha[group] + X * beta);
}
```

Run with `cmdstanr`. Save the full chain including `alpha`, `mu`, `sigma`.

---

## R Package Structure

```
fibr/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── detect_loops.R       # find approximate loops in base space
│   ├── transport_map.R      # estimate H from loop pairs
│   ├── holonomy.R           # main diagnostic function
│   └── plot_holonomy.R      # visualisation
├── src/                     # C++/Rust backend (later)
├── tests/
│   └── testthat/
└── vignettes/
    └── glmm_example.Rmd
```

Main user-facing function:

```r
holonomy_diagnostic(
  chain,           # draws_array or matrix: rows = iterations, cols = parameters
  base_vars,       # character: names of base space parameters, e.g. c("mu", "sigma")
  fiber_vars,      # character: names of fiber parameters, e.g. paste0("alpha[", 1:8, "]")
  epsilon = NULL,  # loop detection tolerance; NULL = auto (based on chain spread)
  n_bootstrap = 200
)
```

Returns an S3 object of class `fibr_holonomy` with print and plot methods.

---

## Dependencies

- `cmdstanr` — for generating chains
- `posterior` — for working with `draws_array` objects
- `Matrix` — sparse linear algebra for large fibers
- `ggplot2` — plotting
- `FNN` or `RANN` — k-nearest neighbours for loop detection

No compiled code in the first milestone — pure R is fine for the diagnostic prototype.

---

## What Comes After

Once the diagnostic is working and validated on the toy model, the roadmap is:

1. **Connection computation**: Extract the off-diagonal Fisher metric blocks $G_{\theta z}$ from the chain or via AD on the log-posterior
2. **Parallel transport correction**: Use the connection to pre-rotate fiber proposals when hyperparameters move
3. **Horizontal leapfrog integrator**: A modified HMC integrator that respects the horizontal/vertical decomposition of the tangent bundle
4. **RMHMC on the total space**: Full implementation with the Fisher metric on $E$, with R bindings via Rust or C++

The diagnostic is the foundation — it validates the theory and gives us a benchmark to measure improvement against at each subsequent stage.

---

## References

- Girolami & Calderhead (2011). Riemann Manifold Langevin and Hamiltonian Monte Carlo. *JRSS-B*.
- Betancourt (2013). A General Metric for Riemannian Manifold HMC. arXiv:1212.4693
- Zhang & Sutton (2014). Semi-Separable HMC for Inference in Bayesian Hierarchical Models. arXiv:1406.3843
- Dhaka et al. (2021). Log-density gradient covariance and automatic metric tensors for Riemann manifold Monte Carlo. arXiv:2211.01746

---

*Generated from a conversation with Claude (claude.ai) on 2026-05-29. The user is an R programmer. The project is at `/Users/bindoffa/antigravity_projects/fibr`.*
