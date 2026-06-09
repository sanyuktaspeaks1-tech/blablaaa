# Chapter 10: Confidence Intervals and Coverage Probability

---

## Coverage Probability of the Normal CI across $\theta$ Values

For $X_1, \ldots, X_n \sim \mathcal{N}(\theta, \sigma^2)$ with $\sigma$ known, sweeps $\theta \in [-1, 1]$ and approximates the coverage probability of

$$\bar{X}_n \pm z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$$

at each point using $M = 1000$ replications, confirming it equals $1 - \alpha$ uniformly.

```r
theta_vals = seq(-1, 1, by = 0.1)
M = 1000        # number of replications
n = 10          # sample size
alpha = 0.05    # confidence level
sigma = 1       # population standard deviation

conf_prob = numeric(length = length(theta_vals))
for (i in 1:length(theta_vals)) {
  theta = theta_vals[i]
  counter = 0
  for (j in 1:M) {
    x = rnorm(n = n, mean = theta, sd = sigma)
    lower_CI = mean(x) - qnorm(1 - alpha/2) * sigma / sqrt(n)
    upper_CI = mean(x) + qnorm(1 - alpha/2) * sigma / sqrt(n)
    if ((lower_CI < theta) && (theta < upper_CI))
      counter = counter + 1
  }
  conf_prob[i] = counter / M   # proportion of containment
}

print(conf_prob)   # confidence coefficient

par(mfrow = c(1, 1))
plot(theta_vals, conf_prob, type = "p", col = "red",
     cex = 1.4, xlab = expression(theta), ylab = "Coverage Probability",
     pch = 19, ylim = c(0.8, 1), cex.lab = 1.4)
```

---

## Coverage Probability of the Normal CI across Sample Sizes $n$

Fixes $\theta_0 = 0$ and varies $n = 1, \ldots, 500$, showing that the coverage probability of $\bar{X}_n \pm z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$ remains flat at $1 - \alpha = 0.95$ regardless of sample size, since the pivot

$$\frac{\bar{X}_n - \theta}{\sigma/\sqrt{n}} \sim \mathcal{N}(0, 1)$$

exactly.

```r
n_vals  = 1:500   # sample sizes
theta_0 = 0       # fixed value of theta
alpha   = 0.05    # size of the test
sigma   = 1       # population standard deviation
M       = 1000    # number of replications

conf_prob = numeric(length = length(n_vals))
for (n in n_vals) {
  counter = 0
  for (i in 1:M) {
    x = rnorm(n = n, mean = theta_0, sd = sigma)
    lower_CI = mean(x) - qnorm(1 - alpha/2) * sigma / sqrt(n)
    upper_CI = mean(x) + qnorm(1 - alpha/2) * sigma / sqrt(n)
    if ((lower_CI < theta_0) & (theta_0 < upper_CI))
      counter = counter + 1
  }
  conf_prob[n] = counter / M
}

plot(n_vals, conf_prob, type = "l", col = "grey",
     xlab = "sample size (n)", cex.lab = 1.5,
     main = expression(P[theta[0]](CI[n])), lwd = 3,
     cex.main = 1.5, ylab = "Confidence Coefficient",
     ylim = c(0.5, 1))
abline(h = 1 - alpha, col = "red", lwd = 3, lty = 2)
```

---

## Coverage Probability for Exponential CI at Fixed $n$

For $X_i \sim \text{Exp}(1/\lambda_0)$ with $\lambda_0 = 3$, constructs the asymptotic CI

$$\bar{X}_n \pm z_{\alpha/2}\frac{\bar{X}_n}{\sqrt{n}}$$

using the estimated SE, and approximates coverage at $n = 10{,}000$ over $M = 1000$ replications. Normal quantiles are justified by the CLT:

$$\frac{\bar{X}_n - \lambda}{\lambda/\sqrt{n}} \xrightarrow{d} \mathcal{N}(0, 1).$$

```r
alpha    = 0.05
n        = 10000   # sample size
lambda_0 = 3       # true value
M        = 1000    # number of replications
counter  = 0       # incremented when CI contains true value

for (i in 1:M) {
  x = rexp(n = n, rate = 1 / lambda_0)
  lower_CI = mean(x) - qnorm(1 - alpha/2) * mean(x) / sqrt(n)
  upper_CI = mean(x) + qnorm(1 - alpha/2) * mean(x) / sqrt(n)
  if ((lower_CI < lambda_0) & (lambda_0 < upper_CI))
    counter = counter + 1
}

conf_prob = counter / M
print(conf_prob)
```

---

## Coverage Probability of Exponential CI Converging to $1 - \alpha$ as $n \to \infty$

Varies $n = 1, \ldots, 500$ for the same asymptotic CI, showing the coverage probability converges to the nominal level $1 - \alpha = 0.95$ as $n$ grows, illustrating the asymptotic confidence coefficient.

```r
n_vals   = 1:500   # varying sample size
lambda_0 = 3       # true value
M        = 1000    # number of replications

conf_prob = numeric(length = length(n_vals))
for (n in n_vals) {
  counter = 0
  for (i in 1:M) {
    x = rexp(n = n, rate = 1 / lambda_0)
    lower_CI = mean(x) - qnorm(1 - alpha/2) * mean(x) / sqrt(n)
    upper_CI = mean(x) + qnorm(1 - alpha/2) * mean(x) / sqrt(n)
    if ((lower_CI < lambda_0) & (lambda_0 < upper_CI))
      counter = counter + 1
  }
  conf_prob[n] = counter / M
}

plot(n_vals, conf_prob, type = "l", col = "grey",
     xlab = "sample size (n)", cex.lab = 1.5,
     main = expression(P[lambda[0]](CI[n])), lwd = 2,
     cex.main = 1.5, ylab = "Confidence Coefficient")
abline(h = 1 - alpha, col = "red", lwd = 3, lty = 2)
```

---

## Coverage Probability of $[a + X_{(n)},\ b + X_{(n)}]$ Depends on $\theta$ — Uniform$(0, \theta)$

For $X_i \sim \text{Uniform}(0, \theta)$, sweeps $\theta \in [0.5, 3]$ and approximates the coverage of $[a + X_{(n)},\ b + X_{(n)}]$ via $M = 5000$ replications, then overlays the exact probability

$$P_\theta = \left(1 - \frac{a}{\theta}\right)^n - \left(1 - \frac{b}{\theta}\right)^n,$$

confirming the coverage depends on the unknown parameter $\theta$.

```r
n          = 10                        # sample size
theta_vals = seq(0.5, 3, by = 0.05)   # discretised parameter space
M          = 5000                      # number of replications
a          = 0.1
b          = 0.2

conf_prob = numeric(length = length(theta_vals))
for (i in 1:length(theta_vals)) {
  theta   = theta_vals[i]
  counter = 0
  for (j in 1:M) {
    x        = runif(n = n, min = 0, max = theta)
    lower_CI = a + max(x)
    upper_CI = b + max(x)
    if ((lower_CI < theta) & (theta < upper_CI))
      counter = counter + 1
  }
  conf_prob[i] = counter / M
}

plot(theta_vals, conf_prob, type = "p", pch = 19,
     col = "red", xlab = expression(theta), cex.lab = 1.5,
     main = bquote(paste("(a,b)", " = (", .(a), ",", .(b), ")")),
     cex.main = 1.4,
     ylab = expression(P[theta](a + X[(n)], b + X[(n)])))

# Overlay exact probability
exact_prob = numeric(length = length(theta_vals))
for (i in 1:length(theta_vals)) {
  exact_prob[i] = (1 - a / theta_vals[i])^n - (1 - b / theta_vals[i])^n
}
lines(theta_vals, exact_prob, type = "l", col = "blue", lwd = 2)
```

---

## Coverage Probability of $[cX_{(n)},\ dX_{(n)}]$ is Constant in $\theta$ — Uniform$(0, \theta)$

For the pivotal CI $[cX_{(n)},\ dX_{(n)}]$ with $1 \leq c < d$, sweeps $\theta \in [0.5, 3]$ and confirms via simulation that the coverage

$$\left(\frac{1}{c}\right)^n - \left(\frac{1}{d}\right)^n$$

is free of $\theta$, yielding a genuine confidence coefficient.

```r
n          = 10                        # sample size
theta_vals = seq(0.5, 3, by = 0.05)   # different theta values
M          = 5000                      # number of replications
c          = 1                         # value of c
d          = 1.3                       # value of d

conf_prob = numeric(length = length(theta_vals))
for (i in 1:length(theta_vals)) {
  theta   = theta_vals[i]
  counter = 0
  for (j in 1:M) {
    x        = runif(n = n, min = 0, max = theta)
    lower_CI = c * max(x)   # lower CI
    upper_CI = d * max(x)   # upper CI
    if ((lower_CI < theta) & (theta < upper_CI))
      counter = counter + 1   # update the counter
  }
  conf_prob[i] = counter / M   # compute the proportion
}

plot(theta_vals, conf_prob, type = "p", pch = 19,
     col = "red", xlab = expression(theta), cex.lab = 1.2,
     main = bquote(paste("(c,d)", " = (", .(c), ",", .(d), ")")),
     cex.main = 1.4, ylab = expression(P[theta](CI[n])),
     ylim = c(0, 1))
```

---

## Confidence Interval via Test-Statistic Inversion — $\mathcal{N}(\mu, \sigma^2)$ with $\sigma$ Known

Visualises how inverting the acceptance region

$$|\bar{x} - \mu_0| \leq z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$$

of a level-$\alpha$ test for $H_0: \mu = \mu_0$ yields the CI $\bar{X}_n \pm z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$. Plots the two boundary lines $\bar{x} = \mu \pm z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$, marks the acceptance region $A(\mu_0)$, and shows the resulting interval $C(\bar{X})$ for a simulated $\bar{X}$.

```r
alpha  = 0.05                    # size of the test
cut_off = 2 * pnorm(alpha / 2)   # one-sample z-test cutoff
sigma  = 2                       # population sd (known)
n      = 5                       # sample size
mu_0   = 2                       # null hypothesis value

plot(0, 0, xlim = c(-2, 7), ylim = c(-1, 6),
     xlab = expression(mu), ylab = expression(bar(x)), cex = 0.2)
abline(h = 0, lwd = 3, col = "grey")
abline(v = 0, lwd = 3, col = "grey")

curve(x - cut_off * sigma / sqrt(n), add = TRUE, col = "red",  lwd = 3, lty = 1)
curve(x + cut_off * sigma / sqrt(n), add = TRUE, col = "red",  lwd = 3, lty = 1)

points(mu_0, 0, pch = 19, col = "magenta")
text(mu_0, -0.5, expression(mu[0]))

segments(mu_0, 0,
         mu_0, mu_0 + cut_off * sigma / sqrt(n),
         lwd = 3, col = "magenta", lty = 2)
segments(mu_0, mu_0 + cut_off * sigma / sqrt(n),
         0,   mu_0 + cut_off * sigma / sqrt(n),
         lwd = 3, col = "magenta", lty = 2)
segments(mu_0, mu_0 - cut_off * sigma / sqrt(n),
         0,   mu_0 - cut_off * sigma / sqrt(n),
         lwd = 3, col = "magenta", lty = 2)

text(-0.5, mu_0, expression(A(mu[0])), col = "magenta")

x     = rnorm(n = n, mean = 5, sd = sigma)
x_bar = mean(x)

text(-0.5, x_bar, expression(bar(X)))
points(0, x_bar, pch = 19, col = "blue")

segments(0, x_bar,
         x_bar + cut_off * sigma / sqrt(n), x_bar,
         col = "blue", lwd = 4, lty = 3)
segments(x_bar - cut_off * sigma / sqrt(n), 0,
         x_bar - cut_off * sigma / sqrt(n), x_bar,
         lwd = 4, col = "blue", lty = 3)
segments(x_bar + cut_off * sigma / sqrt(n), 0,
         x_bar + cut_off * sigma / sqrt(n), x_bar,
         lwd = 4, col = "blue", lty = 3)

text(x_bar, -0.5, expression(C(bar(X))), col = "blue")
```

---

## Visualisation of 100 Repeated 95% CIs for $\mathcal{N}(\mu = 0, 1)$

Simulates 100 independent datasets of size $n = 20$, constructs a 95% CI from each, and plots all intervals — colouring in dark grey those that contain the true $\mu = 0$ and in red those that miss it, illustrating that approximately 5 out of 100 intervals fail to cover the true value.

```r
rep    = 100   # number of replications
sigma2 = 1     # population variance
alpha  = 0.05
n      = 20    # sample size

lower_CI = numeric(rep)   # initialisation
upper_CI = numeric(rep)   # initialisation

for (i in 1:rep) {
  x          = rnorm(n = n, mean = 0, sd = sqrt(sigma2))
  lower_CI[i] = mean(x) - qnorm(1 - alpha/2) * sqrt(sigma2 / n)
  upper_CI[i] = mean(x) + qnorm(1 - alpha/2) * sqrt(sigma2 / n)
}

plot(1:rep, lower_CI, type = "p", pch = 19, col = "grey",
     ylim = c(min(lower_CI), max(upper_CI)),
     ylab = "95% CI", xlab = "replication")
points(1:rep, upper_CI, type = "p", pch = 19, col = "grey")

for (i in 1:rep) {
  if ((lower_CI[i] < 0) * (0 < upper_CI[i]) == 1) {
    segments(i, lower_CI[i], i, upper_CI[i], col = "darkgrey", lwd = 2)
  } else {
    segments(i, lower_CI[i], i, upper_CI[i], col = "red", lwd = 2)
    points(i, lower_CI[i], type = "p", pch = 19, col = "red")
    points(i, upper_CI[i], type = "p", pch = 19, col = "red")
  }
}

abline(h = 0, lwd = 2, lty = 2, col = "blue")
```
