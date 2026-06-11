# Chapter 13
---

**Bootstrap basics – SRSWR and SRSWOR**

Demonstrates the difference between sampling without and with replacement using R's `sample()` function, the foundation of the bootstrap method.

```r
x = c(1,7,5,4,8)
sample(x)                   # SRSWOR
sample(x, replace = TRUE)   # SRSWR
```

---

**Bootstrap standard error of the median**

Generates B = 100 bootstrap samples from x by SRSWR, computes the median for each, and estimates the bootstrap variance and standard error of the sample median.

$$V_{\text{boot}}(T_n) = \frac{1}{B-1}\sum_{b=1}^{B}\left(T_n^{(b)} - \frac{1}{B}\sum_{i=1}^{B} T_n^{(i)}\right)^2$$

```r
B = 100
T_n = numeric(B)
for(i in 1:B){
  T_n[i] = median(sample(x, replace = TRUE))
}
boot_var = (sum((T_n - mean(T_n))^2))/B
print(boot_var)
cat("The estimated bootstrap standard error is ", sqrt(boot_var))
```

---

**Bootstrap CI – step (a): simulate data**

Simulates two independent samples from N(μ₁=6, σ₁²=4) and N(μ₂=4, σ₂²=16) and computes the observed estimate of θ = μ₁ − μ₂.

$$\theta = \mu_1 - \mu_2, \quad \widehat{\theta} = \overline{X}_n - \overline{Y}_m$$

```r
set.seed(123)
mu_1 = 6
x = round(rnorm(n = 20, mean = mu_1, sd = 2),2)
print(x)

set.seed(12)
mu_2 = 4
y = round(rnorm(n = 20, mean = mu_2, sd = 4),2)
print(y)
theta = mu_1 - mu_2
theta_hat = mean(x) - mean(y)
```

---

**Bootstrap CI – step (b): bootstrap sampling distribution of θ̂**

Generates B = 1000 bootstrap samples from both x and y, computes θ̂* = X̄* − Ȳ* for each, plots the bootstrap distribution, and estimates the bootstrap standard error.

```r
B = 1000
boot_theta_hat = numeric(length = B)
for (i in 1:B) {
  x_star = sample(x, size = length(x), replace = TRUE)
  y_star = sample(y, size = length(y), replace = TRUE)
  boot_theta_hat[i] = mean(x_star)- mean(y_star)
}
hist(boot_theta_hat, probability = TRUE,
     main = expression(paste("Bootstrap distribution of ", hat(theta[n]))),
     xlab = expression(theta=mu[1]-mu[2]), cex.lab = 1.4)
se_boot = sqrt(var(boot_theta_hat))
```

---

**Bootstrap CI – step (c): normality-based confidence interval**

Constructs and overlays the (1−α)% normality-based bootstrap confidence interval using the bootstrap standard error.

$$\left(\widehat{\theta} - z_{\alpha/2}\cdot\widehat{\text{SE}}_{\text{boot}},\;\; \widehat{\theta} + z_{\alpha/2}\cdot\widehat{\text{SE}}_{\text{boot}}\right)$$

```r
alpha = 0.05
Normal_CI = c(theta_hat - qnorm(1-alpha/2)*se_boot, theta_hat + qnorm(1-alpha/2)*se_boot)
print(Normal_CI)
abline(v = Normal_CI, col = "red", lwd = 3, lty = 3)
```

---

**Bootstrap CI – step (d): pivotal confidence interval**

Constructs and overlays the pivotal bootstrap confidence interval using the quantiles of the bootstrap distribution.


```r
Pivotal_CI = c(2*theta_hat - quantile(boot_theta_hat, 1-alpha/2), 2*theta_hat - quantile(boot_theta_hat, alpha/2))
print(Pivotal_CI)
abline(v = Pivotal_CI, col = "blue", lwd=3, lty = 4)
```

---

**Bootstrap CI – step (e): percentile confidence interval**

Constructs and overlays the percentile bootstrap confidence interval directly from the empirical quantiles of the bootstrap distribution, and adds a legend comparing all three methods.


```r
Percentile_CI = c(quantile(boot_theta_hat, alpha/2), quantile(boot_theta_hat, 1-alpha/2))
abline(v = Percentile_CI, col = "magenta", lwd = 3, lty =5)
legend("topright", c("Normal", "Pivotal", "Percentile"),
       col = c("red", "blue", "magenta"), lwd = rep(3,3),
       lty = 3:5, bty = "n")
```

---

**Bootstrap regression – step (a): simulate data**

Simulates n = 101 observations from the population regression model y = 0.5 + x + ε with ε ~ N(0, 0.09) and stores them as a data frame.

$$y_i = \beta_0 + \beta_1 x_i + \epsilon_i, \quad \epsilon_i \sim \mathcal{N}(0, \sigma^2)$$

```r
set.seed(123)                                 # for reproducibility
x = seq(0, 1, by = 0.01)
y = 0.5 + 1*x + rnorm(n = length(x), 0, 0.3)
plot(x,y, type = "p", pch = 19, cex = 1.3,col = "red", cex.lab = 1.4)
data = data.frame(x,y)                        # original data
```

---

**Bootstrap regression – step (b): non-parametric bootstrap**

Generates B = 100 bootstrap datasets by resampling rows with replacement (SRSWR), fits a regression line to each, overlays all fitted lines (magenta) and the original fit (blue), and stores the bootstrap estimates of β₀ and β₁.

```r
n = nrow(data)
B = 100  # number of bootstrap replication
beta_0 = numeric(B)              # bootstrap estimate of $\beta_0$
beta_1 = numeric(B)              # bootstrap estimate of $\beta_1$
par(mfrow = c(1,1))
for(i in 1:B){
  ind = sample(1:n, replace = TRUE)
  d = data[ind, ]                 # bootstrap data set D*
  fit = lm( y ~ x, data = d)      # fitting on bootstrap data
  abline(fit, col = "magenta", lty = 2, lwd = 2)
  beta_0[i]=fit$coefficients[1]  # estimate of beta_0 using D*
  beta_1[i]=fit$coefficients[2]  # estimate of beta_1 using D*
}
abline(lm(y ~ x, data = data), col = "blue", lwd = 2)
```

---

**Bootstrap regression – step (c): sampling distributions of β̂₀ and β̂₁**

Plots histograms of the B = 100 bootstrap estimates of β₀ and β₁ and overlays a normal curve with bootstrap mean and standard error on each, confirming that both distributions are well approximated by the normal distribution.

```r
par(mfrow = c(1,2))
hist(beta_0,probability=TRUE,xlab= expression(widehat(beta[0])),
     cex.lab = 1.4, main = "", breaks = 20)
curve(dnorm(x, mean = mean(beta_0),sd = sd(beta_0)), add = TRUE,
      col = "red", lwd = 2)
hist(beta_1, probability = TRUE, xlab = expression(widehat(beta[1])),
     cex.lab = 1.4, main = "", breaks = 20)
curve(dnorm(x, mean = mean(beta_1), sd = sd(beta_1)),add = TRUE,
      col = "red", lwd = 2)
mean(beta_0)             # bootstrap mean
mean(beta_1)             # bootstrap mean
sd(beta_0)               # bootstrap standard error of beta_0_hat
sd(beta_1)               # bootstrap standard error of beta_1_hat
```

---

**Parametric bootstrap for simple linear regression**

Fits the regression model on the original data, then generates B = 100 parametric bootstrap datasets by simulating new responses from N(Ŷᵢ, σ̂²), refits the model on each, and plots the bootstrap sampling distributions of β̂₀ and β̂₁ with normal overlays.

$$Y_i^* \sim \mathcal{N}\!\left(\widehat{\beta}_0 + \widehat{\beta}_1 X_i,\; \widehat{\sigma}^2\right), \quad 1 \leq i \leq n$$

```r
fit = lm(y ~ x, data = data)      # linear regression fit
B = 100                            # number of bootstrap samples
boot_b0 = numeric(B)               # bootstrap estimates of b0
boot_b1 = numeric(B)               # bootstrap estimates of b1
for(i in 1:B){
  mu = coef(fit)[1] + coef(fit)[2]*x       # mean vector
  sd = sqrt((residuals(fit))^2/nrow(data)) # residual sd
  y_star = rnorm(n = nrow(data), mean = mu, sd = sd)
  boot_fit = lm(y_star ~ x)                # fitting bootstrap data
  boot_b0[i] = coef(boot_fit)[1]           # bootstrap estimate b0
  boot_b1[i] = coef(boot_fit)[2]           # bootstrap estimate b1
}
par(mfrow = c(1,2))
hist(boot_b0, probability = TRUE, xlab = expression(widehat(beta[0])),
     cex.lab = 1.4, main = "", breaks = 20)
curve(dnorm(x, mean = mean(boot_b0), sd = sd(boot_b0)),
add = TRUE, col = "red", lwd = 2)
hist(boot_b1, probability = TRUE, xlab = expression(widehat(beta[1])),
     cex.lab = 1.4, main = "", breaks = 20)
curve(dnorm(x, mean = mean(boot_b1), sd = sd(boot_b1)),
add = TRUE, col = "red", lwd = 2)
mean(boot_b0)        # bootstrap mean
mean(boot_b1)        # bootstrap mean
sd(boot_b0)          # bootstrap standard error of beta_0_hat
sd(boot_b1)          # bootstrap standard error of beta_1_hat
```

---

**Bootstrap mean convergence as B → ∞**

Simulates a fixed sample of n = 5 from N(0,1) and tracks the bootstrap double mean X̄*_{n,B} as B increases from 1 to 500, showing convergence to the sample mean X̄_n rather than the population mean μ.

$$\overline{\overline{X}}_{n,B}^* = \frac{1}{B}\sum_{b=1}^{B}\overline{X}_n^{*(b)} \xrightarrow{B\to\infty} \overline{X}_n \xrightarrow{n\to\infty} \mu$$

```r
set.seed(12)
mu = 0
sigma = 1
x = rnorm(n = 5, mean = mu, sd = sigma)
x_bar = mean(x)
B_vals = 1:500
boot_means = numeric(length = length(B_vals))
for(B in B_vals){
  boot_mean = numeric(B)
  for(i in 1:B){
    x_star = sample(x, replace = TRUE)
    boot_mean[i] = mean(x_star)
  }
  boot_means[B] = mean(boot_mean)
}
par(mfrow = c(1,1))
plot(B_vals, boot_means, type = "l",
     col = "darkgrey", lwd = 2, xlab = "B",
     cex.lab = 1.2, ylim = c(-1, 0.5), cex.main = 1.3)
abline(h = x_bar, col = "blue", lty = 2, lwd = 4)
abline(h = mu, col = "magenta", lwd = 4, lty = 3)
legend("topleft", legend = c(expression(mu), expression(bar(X[n])),
                             expression(bar(bar(X[n*B])))),
       col = c("magenta", "blue", "grey"), lty = c(3,2,1),
       lwd = c(4,4,3), bty = "n", cex = 1.4)
```

---

**Bootstrap for logistic growth – step I: simulate data**

Simulates population size from the logistic growth model with x₀ = 5, r = 0.8, K = 20, σ = 2 as the dataset for bootstrapping.

```r
LogisticModel = function(t, x0, r, K){
  K/(1+(K/x0 - 1)*exp(-r*t))
}
control = nls.control(maxiter = 500)
x0 = 5; K = 20; r = 0.8
t = seq(0, 15, by = 0.3)
x = LogisticModel(t, x0 = 5, r = 0.8, K = 20) + rnorm(length(t), 0, 2)
data = data.frame(t, x)
```

---

**Bootstrap for logistic growth – step II: bootstrap NLS fits**

Generates B = 500 bootstrap datasets by resampling rows with replacement and refits the logistic model via `nls()` on each, storing the bootstrap estimates of x₀, K, and r_m.

```r
B = 500
boot.x0 = numeric(B)              # bootstrap estimate of x0
boot.r = numeric(B)               # bootstrap estimate of r_m
boot.K = numeric(B)               # bootstrap estimate of K
for(i in 1:B){
  brow = sample(1:nrow(data), replace = T)
  newdata = data[brow,]
  fit = nls(x~LogisticModel(t, x0, r, K), control=control,
            data = newdata, start = list(x0 = 5, K = 20, r = 0.8))
  boot.x0[i] = coef(fit)[1]
  boot.K[i] = coef(fit)[2]
  boot.r[i] = coef(fit)[3]
}
```

---

**Bootstrap for logistic growth – step III: sampling distributions**

Plots the bootstrap sampling distributions of x̂₀, K̂, and r̂_m as histograms side by side based on B = 500 bootstrap replications.

```r
par(mfrow=c(1,3))
hist(boot.x0, probability = T, xlab = expression(widehat(x[0])),
     main = paste("B = ", B), breaks = 20, cex.lab = 1.4)
hist(boot.K, probability = T, xlab = expression(widehat(K)),
     main = paste("B = ", B), breaks = 20, cex.lab = 1.4)
hist(boot.r, probability = T, xlab = expression(widehat(r[m])),
     main = paste("B = ", B), breaks = 20, cex.lab  = 1.4)
```

---

**Bootstrap for logistic growth – step IV: pairwise scatterplots**

Plots all pairwise scatterplots of x̂₀, K̂, and r̂_m to reveal the nonlinear dependencies and correlations between the bootstrap parameter estimates, particularly the strong correlation between x̂₀ and r̂_m.

```r
par(mfrow = c(2,3))
plot(boot.x0, boot.K, type = "p", col = "red",xlab = expression(widehat(x[0])),ylab = expression(widehat(K)))
plot(boot.x0, boot.r, type = "p", col = "red",xlab = expression(widehat(x[0])),ylab = expression(widehat(r[m])))
plot(boot.K, boot.x0, type = "p", col = "red",xlab = expression(widehat(K)), ylab = expression(widehat(x[0])))
plot(boot.K, boot.r, type = "p", col="red", xlab = expression(widehat(K)), ylab = expression(widehat(r[m])))
plot(boot.r, boot.x0, type = "p", col="red",xlab = expression(widehat(r[m])), ylab = expression(widehat(x[0])))
plot(boot.r, boot.K, type = "p", col="red", xlab = expression(widehat(r[m])),ylab = expression(widehat(K)))
```

---

**Teacher's corner – bootstrap sampling distribution of the sample mean**

Generates B = 1000 bootstrap samples from x = {1, 2, 3, 4, 5} and plots the histogram of bootstrap sample means, demonstrating to students that the bootstrap distribution appears approximately normal.

```r
x = 1:5                                  # values {1,2,3,4,5}
B = 1000                                 # number of bootstrap samples
boot_means = numeric(B)
for(i in 1:B){
  boot_means[i] = mean(sample(x, size = length(x),replace =  TRUE))
}
hist(boot_means, probability = TRUE,
     xlab = expression(bar(X[n])^B), main = "B = 1000")
```
