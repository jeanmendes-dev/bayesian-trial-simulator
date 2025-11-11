# Bayesian Clinical Trial Simulator

[![R-CMD-check](https://github.com/rforbio/bayesian-trial-simulator/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rforbio/bayesian-trial-simulator/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> A lightweight, transparent R framework for simulating **Bayesian adaptive clinical trials** with sequential posterior updating and early stopping rules — designed for statisticians, clinical teams, and biotech decision-makers.

🎯 **Use cases**:  
- Design exploration for Phase II oncology trials (binary endpoints: ORR, CR)  
- Operating characteristics assessment (power, type I error, sample size distribution)  
- Regulatory submission prep (FDA/EMA encourage Bayesian adaptive designs)  
- Academic teaching & industry training

---

## 🧠 Why Bayesian Adaptive Designs?

Traditional fixed-sample designs often lead to:
- ❌ Over-enrollment when superiority is clear early  
- ❌ Delayed stopping when futility is evident  
- ❌ Inflexibility in rare disease settings (small N, high uncertainty)

Bayesian adaptive trials offer:
- ✅ Real-time evidence updating via Bayes’ theorem  
- ✅ Ethical efficiency: stop early for efficacy/futility  
- ✅ Seamless incorporation of historical/prior knowledge  
- ✅ Direct probability statements clinicians understand:  
  > *“There is a 92% probability the new therapy outperforms control.”*

---

## 🛠️ Core Features

| Component | Description |
|---------|-------------|
| **Sequential Simulation** | Simulate trials coorte-by-coorte (e.g., n=6, 12, 18…) |
| **Conjugate Updating** | Beta-binomial (binary), Normal-Normal (continuous) |
| **Decision Rules** | Customizable thresholds for:<br> • Efficacy: `P(θₜ > θ꜀ + δ \| data) > γ₁`<br> • Futility: `P(θₜ > θ꜀ \| data) < γ₂` |
| **Visualization** | Evolution of posterior, `P(treatment better)`, final densities |
| **Extensible** | Add models (Weibull, hierarchical), priors (empirical, weakly informative), allocation (Thompson sampling) |

---

## ▶️ Quick Example

Simulate a Phase II oncology trial (binary endpoint: objective response rate):

```r
source("https://raw.githubusercontent.com/rforbio/bayesian-trial-simulator/main/bayesian_trial_simulator.R")

trial <- simulate_bayesian_trial(
  n_total = 80,
  cohort_size = 10,
  true_rate_control = 0.30,   # Standard of care
  true_rate_treatment = 0.50, # Novel therapy
  prior_control = c(1, 1),    # Uniform prior
  prior_treatment = c(1, 1),
  efficacy_threshold = 0.90,
  futility_threshold = 0.10,
  seed = 2025
)

print(trial)
plot_trial_evolution(trial)
