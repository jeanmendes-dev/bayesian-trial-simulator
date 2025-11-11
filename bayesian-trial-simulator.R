#if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
#if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
#if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
#if (!requireNamespace("gganimate", quietly = TRUE)) install.packages("gganimate")

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# Função principal: simular ensaio clínico adaptativo binário
# ------------------------------------------------------------
simulate_bayesian_trial <- function(
    n_total = 80,
    cohort_size = 10,
    true_rate_control = 0.30,   # taxa real de resposta no braço controle (ex: padrão)
    true_rate_treatment = 0.50, # taxa real no braço novo (ex: droga experimental)
    prior_control = c(alpha = 1, beta = 1),   # Beta(1,1) = uniforme
    prior_treatment = c(alpha = 1, beta = 1),
    efficacy_threshold = 0.90,  # P(treat > control | data) > 90% → parar por eficácia
    futility_threshold = 0.10,  # P(treat > control | data) < 10% → parar por futilidade
    seed = 123
) {
  set.seed(seed)
  
  n_cohorts <- ceiling(n_total / cohort_size)
  
  # Armazenar histórico
  history <- data.frame(
    cohort = integer(),
    n_control = integer(),
    n_treatment = integer(),
    y_control = integer(),
    y_treatment = integer(),
    post_alpha_c = numeric(),
    post_beta_c = numeric(),
    post_alpha_t = numeric(),
    post_beta_t = numeric(),
    prob_treat_better = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  
  # Inicializar contadores
  n_c <- n_t <- 0
  y_c <- y_t <- 0
  post_alpha_c <- prior_control[1]
  post_beta_c  <- prior_control[2]
  post_alpha_t <- prior_treatment[1]
  post_beta_t  <- prior_treatment[2]
  
  for (k in 1:n_cohorts) {
    # Alocar metade do coorte para cada braço (ou ajuste para randomização adaptativa depois)
    n_cohort_half <- floor(cohort_size / 2)
    n_c_new <- n_cohort_half
    n_t_new <- cohort_size - n_cohort_half
    
    # Gerar dados simulados (verdade desconhecida para análise)
    y_c_new <- rbinom(1, n_c_new, true_rate_control)
    y_t_new <- rbinom(1, n_t_new, true_rate_treatment)
    
    # Atualizar contagens
    n_c <- n_c + n_c_new
    n_t <- n_t + n_t_new
    y_c <- y_c + y_c_new
    y_t <- y_t + y_t_new
    
    # Atualizar posterior (Beta-Binomial conjugado)
    post_alpha_c <- prior_control[1] + y_c
    post_beta_c  <- prior_control[2] + n_c - y_c
    post_alpha_t <- prior_treatment[1] + y_t
    post_beta_t  <- prior_treatment[2] + n_t - y_t
    
    # Estimar P(theta_t > theta_c | dados) via Monte Carlo (rápido e robusto)
    nsim <- 10000
    theta_c_sim <- rbeta(nsim, post_alpha_c, post_beta_c)
    theta_t_sim <- rbeta(nsim, post_alpha_t, post_beta_t)
    prob_treat_better <- mean(theta_t_sim > theta_c_sim)
    
    # Regra de decisão
    decision <- "Continue"
    if (prob_treat_better > efficacy_threshold) {
      decision <- "Stop: Efficacy"
    } else if (prob_treat_better < futility_threshold) {
      decision <- "Stop: Futility"
    }
    
    # Registrar no histórico
    new_row <- data.frame(
      cohort = k,
      n_control = n_c,
      n_treatment = n_t,
      y_control = y_c,
      y_treatment = y_t,
      post_alpha_c = post_alpha_c,
      post_beta_c = post_beta_c,
      post_alpha_t = post_alpha_t,
      post_beta_t = post_beta_t,
      prob_treat_better = prob_treat_better,
      decision = decision
    )
    history <- bind_rows(history, new_row)
    
    # Parada antecipada?
    if (decision != "Continue") break
  }
  
  # Adicionar coluna de tamanho total até o momento
  history$total_n <- history$n_control + history$n_treatment
  
  # Retornar
  structure(
    list(
      history = history,
      parameters = list(
        n_total = n_total,
        cohort_size = cohort_size,
        true_rate_control = true_rate_control,
        true_rate_treatment = true_rate_treatment,
        prior_control = prior_control,
        prior_treatment = prior_treatment,
        efficacy_threshold = efficacy_threshold,
        futility_threshold = futility_threshold
      ),
      stopped_early = any(history$decision != "Continue"),
      final_decision = tail(history$decision, 1)
    ),
    class = "bayes_trial"
  )
}

# ------------------------------------------------------------
# Método print para objetos bayes_trial
# ------------------------------------------------------------
print.bayes_trial <- function(x, ...) {
  cat("✅ Bayesian Adaptive Trial Simulation\n")
  cat("====================================\n")
  p <- x$parameters
  cat(sprintf("- True rates: Control = %.0f%% | Treatment = %.0f%%\n", 
              100 * p$true_rate_control, 100 * p$true_rate_treatment))
  cat(sprintf("- Priors: Control = Beta(%.0f, %.0f) | Treatment = Beta(%.0f, %.0f)\n", 
              p$prior_control[1], p$prior_control[2],
              p$prior_treatment[1], p$prior_treatment[2]))
  cat(sprintf("- Stopping rules: Efficacy > %.0f%% | Futility < %.0f%%\n", 
              100 * p$efficacy_threshold, 100 * p$futility_threshold))
  cat("\n📈 Trial Progress:\n")
  print(x$history[, c("cohort", "total_n", "y_control", "y_treatment", 
                      "prob_treat_better", "decision")], row.names = FALSE)
  cat("\n🔍 Final Decision:", x$final_decision, "\n")
  if (x$stopped_early) cat("⚠️  Trial stopped early at n =", tail(x$history$total_n, 1), "\n")
}

# ------------------------------------------------------------
# Função para plotar evolução da probabilidade de superioridade
# ------------------------------------------------------------
plot_trial_evolution <- function(trial_obj, animate = FALSE) {
  df <- trial_obj$history
  
  p <- ggplot(df, aes(x = total_n, y = prob_treat_better)) +
    geom_line(color = "#2C3E50", size = 1.1) +
    geom_point(aes(color = decision), size = 2.5) +
    geom_hline(yintercept = trial_obj$parameters$efficacy_threshold, 
               linetype = "dashed", color = "#27AE60", size = 0.8) +
    geom_hline(yintercept = trial_obj$parameters$futility_threshold, 
               linetype = "dashed", color = "#E74C3C", size = 0.8) +
    scale_color_manual(values = c(
      "Continue" = "#3498DB",
      "Stop: Efficacy" = "#27AE60",
      "Stop: Futility" = "#E74C3C"
    )) +
    labs(
      title = "Adaptive Bayesian Trial: Evolution of P(Treatment > Control)",
      subtitle = "Decision thresholds: Efficacy (green) / Futility (red)",
      x = "Total Sample Size (n)",
      y = "P(θₜ > θ꜀ | Data)",
      color = "Decision"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  if (animate && requireNamespace("gganimate", quietly = TRUE)) {
    p <- p + 
      transition_reveal(total_n) +
      shadow_wake(wake_length = 0.05)
  }
  
  print(p)
  if (animate) {
    cat("\nℹ️  To save animation: anim_save('trial_evolution.gif', last_plot())\n")
  }
}

# ------------------------------------------------------------
# ✅ EXEMPLO DE USO — Oncologia (Fase II, ORR)
# ------------------------------------------------------------
# Simular ensaio: droga nova vs. quimioterapia padrão
# - Padrão: 30% de resposta objetiva (ORR)
# - Nova droga: esperamos 50% (mais realista: 40–55%)
# - Coortes de 10 pacientes (5 por braço)
# - Parar se P(superioridade) > 90% (eficácia) ou < 10% (futilidade)

trial <- simulate_bayesian_trial(
  n_total = 80,
  cohort_size = 10,
  true_rate_control = 0.30,
  true_rate_treatment = 0.50,
  prior_control = c(1, 1),
  prior_treatment = c(1, 1),
  efficacy_threshold = 0.90,
  futility_threshold = 0.10,
  seed = 2025
)

# Mostrar resumo
trial

# Plotar evolução
plot_trial_evolution(trial)

# (Opcional) Animação
# plot_trial_evolution(trial, animate = TRUE)

# ------------------------------------------------------------
# 📊 Análise final: distribuições posteriores
# ------------------------------------------------------------
final <- tail(trial$history, 1)

# Gerar curvas de densidade posterior
theta_grid <- seq(0, 1, length.out = 500)
post_c <- dbeta(theta_grid, final$post_alpha_c, final$post_beta_c)
post_t <- dbeta(theta_grid, final$post_alpha_t, final$post_beta_t)

dens_df <- data.frame(
  theta = rep(theta_grid, 2),
  density = c(post_c, post_t),
  arm = rep(c("Control", "Treatment"), each = length(theta_grid))
)

p <- ggplot(dens_df, aes(x = theta, y = density, fill = arm, color = arm)) +
  geom_area(alpha = 0.3, position = "identity") +
  geom_line(size = 1) +
  labs(
    title = "Final Posterior Distributions",
    subtitle = sprintf("Control: Beta(%d, %d) | Treatment: Beta(%d, %d)",
                       final$post_alpha_c, final$post_beta_c,
                       final$post_alpha_t, final$post_beta_t),
    x = expression(theta ~ "(Response Rate)"),
    y = "Density"
  ) +
  scale_fill_manual(values = c("#E74C3C", "#27AE60")) +
  scale_color_manual(values = c("#C0392B", "#229954")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p)  # 👈 Isso força a renderização no painel Plots!

# Salvar gráfico em alta resolução
ggsave("final_posterior_distributions.png", 
       p, 
       width = 10, 
       height = 6, 
       dpi = 300,
       bg = "white")

cat("\n✅ Gráfico salvo como 'final_posterior_distributions.png'\n")
# ------------------------------------------------------------
# - Expandir para modelos contínuos (ex: PFS com normal)
# - Adicionar randomização adaptativa (Thompson Sampling)
# - Criar pacote R (`usethis::create_package()`)
# - Construir app Shiny (ui/server prontos em próximo passo)
# ------------------------------------------------------------