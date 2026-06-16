# verify_flat_connection.R
# -----------------------------------------------------------------------------
# Claim under test: the Ehresmann connection A = -G_FF^{-1} G_BF for the centred
# logistic GLMM has IDENTICALLY ZERO curvature (it is flat), so its true
# geometric holonomy is trivial. The non-zero F_j = -2/(sigma^5 G_FF^2) reported
# in the paper is the curvature of the LINEARISED connection (G_FF frozen at the
# loop centre / alpha held fixed), which is what synthetic_holonomy_loop()
# integrates.
#
# Two independent checks:
#   (1) Symbolic-style algebra (documented below).
#   (2) Numerical: integrate the transport ODE around a closed loop in (mu,sigma)
#       two ways -- (a) freezing G at the centre (linearised), (b) updating G
#       along the fiber (true connection) -- and compare the net displacement.
#
# Expected result: (a) grows with loop area; (b) is ~0 to integrator precision.
# -----------------------------------------------------------------------------
#
# ALGEBRA (full Ehresmann curvature of a nonlinear connection on a 1-D fiber):
#   For horizontal lift  dalpha = A_mu dmu + A_sigma dsigma, the curvature is
#     F = d_mu A_sigma - d_sigma A_mu  +  A_mu d_alpha A_sigma - A_sigma d_alpha A_mu
#   The first pair is the paper's F_j = -2/(sigma^5 G^2).
#   The second pair (the vertical / fiber-derivative terms omitted in Remark 4)
#   evaluates to +2/(sigma^5 G^2): the dS/dalpha terms cancel exactly, leaving
#   exactly the negative of the first pair. Hence F == 0. (Derivation in the
#   handoff; reproduce with any CAS to confirm.)

set.seed(1)

# One group: n observations with fixed x*beta offsets, prior N(mu, sigma^2)
n  <- 3L                      # sparse / prior-dominated
xb <- rnorm(n, 0, 1)          # x_i^T beta offsets (held fixed)

S      <- function(a) { p <- plogis(a + xb); sum(p * (1 - p)) }          # likelihood Fisher info
G_FF   <- function(a, s) 1 / s^2 + S(a)
A_mu   <- function(a, mu, s) 1 / (s^2 * G_FF(a, s))
A_sig  <- function(a, mu, s) 2 * (a - mu) / (s^3 * G_FF(a, s))

# Holonomy = net alpha displacement after one anticlockwise circular loop in
# (mu, sigma) of radius r, integrated with RK4. `frozen = TRUE` reproduces
# synthetic_holonomy_loop() (G fixed at centre); FALSE is the true connection.
holonomy <- function(mu0, sigma0, alpha0, r, frozen, nsteps = 20000L) {
  G0  <- G_FF(alpha0, sigma0)
  f   <- function(tk, ak) {
    m  <- mu0 + r * cos(tk); s <- sigma0 + r * sin(tk)
    dm <- -r * sin(tk);      ds <- r * cos(tk)
    if (frozen) {
      am <- 1 / (s^2 * G0); asg <- 2 * (ak - m) / (s^3 * G0)
    } else {
      am <- A_mu(ak, m, s); asg <- A_sig(ak, m, s)
    }
    am * dm + asg * ds
  }
  th <- seq(0, 2 * pi, length.out = nsteps + 1L)
  a  <- alpha0
  for (k in seq_len(nsteps)) {
    dt <- th[k + 1L] - th[k]; t <- th[k]
    k1 <- f(t, a); k2 <- f(t + dt / 2, a + dt / 2 * k1)
    k3 <- f(t + dt / 2, a + dt / 2 * k2); k4 <- f(t + dt, a + dt * k3)
    a  <- a + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  a - alpha0
}

mu0 <- 0; sigma0 <- 0.5; alpha0 <- 0.7
cat(sprintf("%-6s %22s %22s\n", "r", "frozen-G (linearised)", "full A (true)"))
for (r in c(0.05, 0.1, 0.2)) {
  hf <- holonomy(mu0, sigma0, alpha0, r, frozen = TRUE)
  hv <- holonomy(mu0, sigma0, alpha0, r, frozen = FALSE)
  cat(sprintf("%-6.2f %22.6e %22.6e\n", r, hf, hv))
}
# Python reference output (RK4, nsteps=20000):
#   r=0.05  frozen=+2.54e-02   full=-1.0e-14
#   r=0.10  frozen=+1.27e-01   full=+4.2e-15
#   r=0.20  frozen=+1.49e+00   full=+1.7e-15
