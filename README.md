# Chapter 14 – Central Limit Theorem and Delta Method Simulations

This repository contains R code that visually demonstrates key asymptotic results from Chapter 14: the Central Limit Theorem (CLT) applied to Poisson sample means, several variants of the Delta method, Wald test statistics, and the multivariate Delta method for a coefficient of variation. Each script simulates sampling distributions, overlays the corresponding normal approximation, and shows how the approximation improves with increasing sample size.

---

## 1. CLT for a Poisson Sample Mean


Simulates **M = 500** replications of the sample mean  ${\overline{X}_n}$ drawn from a $\text{Poisson}(\lambda=3)$ distribution. The sample sizes considered are $n \in \{2, 3, 5, 10, 15, 25\}$. For each $n$, a histogram of the simulated means is plotted and overlaid with the theoretical $\mathcal{N}\bigl(\lambda, \lambda/n\bigr)$ density, illustrating the convergence in distribution stated by the CLT.

$$
\bar{X}_n \sim \mathcal{N}\left(\lambda, \frac{\lambda}{n}\right) \quad \text{for large } n
$$


```r
lambda = 3 # true value of lambda 
par(mfrow = c(2,3))
n_vals = c(2,3,5,10,15,25) # varying sample size
M = 500 # number of replications
for(n in n_vals){ # loop on sample size
  sample_means = numeric(length = M)
  for(i in 1:M){ # loop on number of replications
    x = rpois(n = n, lambda = lambda)
    sample_means[i] = mean(x)
  }
  hist(sample_means, probability = TRUE, main = paste("n = ", n),
       xlab = expression(bar(X[n])), breaks = 20,
       cex.lab = 1.4, cex.main = 1.4)
  curve(dnorm(x, mean = lambda, sd = sqrt(lambda/n)),
        add = TRUE, col = "red", lwd = 2)
} 
```

---

## 2. Delta Method – Sampling Distribution of $\(1/\overline{X}_n\)$

This section explores the sampling distribution of the reciprocal of the sample mean, $\overline{Y}_n = 1/\overline{X}_n$, again using $\text{Poisson}(\lambda=3)$ data. For the same sample sizes $n \in \{2,3,5,10,15,25\}$ and 500 replications, the Delta method provides the asymptotic normal approximation

$$
Y_n = \frac{1}{\overline{X}_n} \sim \mathcal{N}\\left(\frac{1}{\lambda}, \frac{1}{\lambda^3 n}\right) \quad \text{for large } n.
$$

The red curve is the approximating normal density, and the blue dot marks the true value $1/\lambda$.

```r
lambda = 3 # true value of lambda
par(mfrow = c(2,3))
n_vals = c(2,3,5,10,15,25) # varying sample size
M = 500 # number of replications
for(n in n_vals){
  y_n = numeric(length = M)
  for(i in 1:M){
    x = rpois(n = n, lambda = lambda)
    y_n[i] = 1/mean(x)
  }
  hist(y_n, probability = TRUE, main = paste("n = ", n),
       xlab = expression(bar(Y[n])), breaks = 20,
       cex.lab = 1.4, cex.main = 1.4)
  curve(dnorm(x, mean = 1/lambda, sd = sqrt(1/(lambda^3*n))),
        add = TRUE, col = "red", lwd = 2)
  points(1/lambda,0, pch = 19, col = "blue", cex = 1.5)
}
```

---

## 3. Delta Method – Sampling Distribution of $\(\psi(\bar{X}_n)\)$

Here we apply the Delta method to the function $\psi(\bar{X}_n) = 1 - (1+\bar{X}_n)e^{-\bar{X}_n}$, which estimates $\psi(\lambda) = P(X \ge 1)$ for a Poisson distribution with rate $\lambda=2$. The simulation uses sample sizes $n \in \{4, 10, 20, 50, 100, 500\}$ and 1000 replications each. The asymptotic variance is approximated by

$$
\text{Var}(\psi(\bar{X}_n)) \approx \frac{\lambda^3 e^{-2\lambda}}{n}.
$$

The true parameter $\psi(\lambda)$ is shown as a red dot, while the red curve represents the normal approximation.

```r
par(mfrow = c(2,3))
n_vals = c(4, 10, 20, 50, 100, 500)
rep = 1000
lambda = 2
psi = 1-(1+lambda)*exp(-lambda)
for(n in n_vals){
  psi_vals = numeric(length = rep)
  for(i in 1:rep){
    x = rpois(n = n, lambda = lambda)
    psi_vals[i] = 1-(1+mean(x))*exp(-mean(x))
  }
  hist(psi_vals, probability = TRUE, cex.lab=1.5,
       xlab = expression(psi), main = paste("n = ", n))
  points(psi, 0, pch = 19, col = "red", cex = 2)
  curve(dnorm(x,mean = psi,sd=sqrt(lambda^3*exp(-2*lambda)/n)), add = TRUE, col = "red", lwd =2)
}
```

---

## 4. CLT and Delta Method – Wald Test Statistics $\(W_\lambda\)$ and $\(W_\psi\)$

To verify that Wald statistics follow an approximate standard normal distribution under the null, we simulate 1000 replications for sample sizes $n \in \{5, 10, 25\}$ and true parameter $\lambda=4$. Two Wald statistics are constructed:

- $W_\lambda = \frac{\bar{X}_n - \lambda_0}{\widehat{\text{SE}}(\bar{X}_n)}$ using the CLT,
- $W_\psi = \frac{\hat{\psi}_n - \psi_0}{\widehat{\text{SE}}(\hat{\psi}_n)}$ with $\psi = \lambda^2$ using the Delta method, where $\widehat{\text{SE}}(\hat{\psi}_n) = 2\sqrt{\bar{X}_n^3/n}$.
Histograms of both statistics are overlaid with the $\(\mathcal{N}(0,1)\)$ density.

```r
par(mfrow = c(2,3))
n_vals = c(5, 10, 25) # varied sample size
lambda = 4 # true parameter
psi = lambda^2 # true value of psi

for(n in n_vals){
  rep = 1000 # number of replications
  lambda_hat = numeric(rep) # estimate of lambda
  se_lambda_hat = numeric(rep) # estimated SE(lambda_hat)
  psi_hat = numeric(rep) # estimated SE(psi_hat)
  se_psi_hat = numeric(rep)
  w_lambda = numeric(length = rep) # Wald statistics (lambda)
  w_psi = numeric(length = rep) # Wald statistics (psi)
  
  for(i in 1:rep){
    x = rpois(n = n, lambda = lambda)
    lambda_hat[i] = mean(x)
    se_lambda_hat[i] = sqrt(lambda_hat[i]/n)
    psi_hat[i] = lambda_hat[i]^2
    se_psi_hat[i] = 2*sqrt(lambda_hat[i]^3/n)
    
    w_lambda[i] = (lambda_hat[i] - lambda)/se_lambda_hat[i]
    w_psi[i] = (psi_hat[i]-psi)/se_psi_hat[i]
  }
  
  hist(w_lambda, probability = TRUE, main = paste("n = ", n),
       xlab = expression(w[lambda]), breaks = 30,
       cex.lab = 1.5, cex.main = 1.3)
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  
  hist(w_psi, probability = TRUE, main = paste("n = ", n),
       xlab = expression(w[psi]), breaks = 30,
       cex.lab = 1.5, cex.main = 1.3)
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
}
```

---

## 5. Multivariate Delta Method – Sampling Distribution of the Coefficient of Variation

The coefficient of variation $\text{CV} = S_n / \bar{X}_n$ is estimated from i.i.d. $\mathcal{N}(\mu=3, \sigma^2=4)$ data. For sample sizes $n \in \{5, 10, 20, 50, 100, 250\}$ and 1000 replications, the multivariate Delta method yields the approximate variance

$$
\text{Var}\\left(\frac{\hat{\sigma}}{\hat{\mu}}\right) \approx \frac{\sigma^4}{\mu^4} \cdot \frac{\sigma^2}{n} + \frac{1}{\mu^2}\,\text{Var}(\hat{\sigma}),
$$
where
$$
\text{Var}(\hat{\sigma}) = \sigma^2\left[1 - \left(\sqrt{\frac{2}{n-1}}\cdot\frac{\Gamma(n/2)}{\Gamma((n-1)/2)}\right)^2\right].
$$

The normal approximation (red curve) is superimposed on the histogram of the simulated CV values.

```r
n_vals = c(5,10,20,50,100,250) # sample size
M = 1000 # number of replications
mu = 3 # true mean
sigma = 2 # true sd
par(mfrow = c(2,3))
for(n in n_vals){
  sample_cv = numeric(length = M)
  for(i in 1:M){
    x = rnorm(n = n, mean = mu, sd = sigma)
    sample_cv[i] = sd(x)/mean(x) # coefficient of variation
  }
  hist(sample_cv, probability = TRUE, 
       xlab = expression(widehat(CV)), main = paste("n = ",n))
  delta_var = (sigma/mu)^4/n + ((1/mu)^2)*sigma^2*(1-
(sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2))^2)
  curve(dnorm(x, mean = sigma/mu, sd = sqrt(delta_var)),
        add = TRUE, col = "red", lwd = 2)
}
```
