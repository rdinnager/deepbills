# Bayesian Interpretation of Riemannian Manifold Evolutionary Models

## A Mathematical Framework for Modeling Bird Beak Evolution on VAE Latent Manifolds

*Russell Dinnage, with mathematical notes prepared March 2026*

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Background: Brownian Motion in Phylogenetics and on Manifolds](#2-background-brownian-motion-in-phylogenetics-and-on-manifolds)
3. [The VAE's Riemannian Geometry](#3-the-vaes-riemannian-geometry)
4. [Path Energy vs Path Length](#4-path-energy-vs-path-length)
5. [The Onsager-Machlup Functional](#5-the-onsager-machlup-functional)
6. [MAP Estimation: From Loss Function to Posterior](#6-map-estimation-from-loss-function-to-posterior)
7. [What Changed from v2 to v3](#7-what-changed-from-v2-to-v3)
8. [The Rho Annealing Schedule as Posterior Tempering](#8-the-rho-annealing-schedule-as-posterior-tempering)
9. [Gradient Accumulation for Optimization](#9-gradient-accumulation-for-optimization)
10. [Summary and Open Questions](#10-summary-and-open-questions)
11. [References](#11-references)

---

## 1. Introduction

### The Problem

In Dinnage & Kleineberg (2025), we trained a DeepSDF model on 3D beak scans from 2,020 bird species, learning 64-dimensional latent representations that capture ecologically meaningful variation. A conditional VAE (CVAE) was then trained on these representations, mapping beak shape codes and trophic niche information into a joint latent space.

The next question is: **how does beak morphology evolve across the bird phylogeny?** The standard approach would be to fit a Brownian motion (BM) model to the latent codes, treating the 16-dimensional active latent space as flat Euclidean space. But a VAE's latent space is *not* flat. The posterior distribution of the VAE defines a natural Riemannian geometry where some regions are densely populated by real beak shapes and others are "dead zones" that map to implausible morphologies. Treating the space as Euclidean ignores this structure entirely.

Our approach is to model evolution as Brownian motion *on the Riemannian manifold* induced by the VAE. Evolutionary paths along phylogenetic branches are constrained to stay in regions of the latent space that correspond to biologically plausible beak shapes, and the "cost" of traversing different regions varies according to the local data density.

### Why a Bayesian Interpretation Matters

The model is currently fitted by minimizing a composite loss function with several terms weighted by hand-chosen constants. This works well in practice, but it raises questions:

- What *exactly* are our priors?
- How strong or weak are they?
- Can we interpret the loss weights in terms of meaningful statistical quantities?
- Could we, in principle, do full posterior inference (not just point estimation)?

This document shows that, with one key modification to the loss function (using path **energy** instead of path **length**), the entire optimization can be understood as **Maximum A Posteriori (MAP) estimation** under a well-defined Bayesian model. Each loss term corresponds to either a likelihood or a prior, and each weight maps to an interpretable precision (inverse-variance) parameter.

---

## 2. Background: Brownian Motion in Phylogenetics and on Manifolds

### 2.1 Brownian Motion in Phylogenetics

The use of Brownian motion to model trait evolution on phylogenies was introduced by Felsenstein (1985), one of the most cited papers in evolutionary biology. The idea is simple but powerful:

**In plain language:** Imagine a trait (say, beak length) evolving along a branch of the phylogenetic tree. At each instant, the trait value gets a tiny random "kick" — a small random change drawn from a normal distribution. Over time, these kicks accumulate, so the total change over a branch of length *T* is normally distributed with variance proportional to *T*.

**Mathematically:** A continuous trait *z(t)* evolves according to the stochastic differential equation (SDE):

$$dz = \sigma \, dW_t$$

where *W_t* is a standard Wiener process (the mathematical formalization of "random kicks") and *sigma* is the **evolutionary rate** parameter controlling how fast the trait changes. Over a branch of duration *T*, the displacement is:

$$z(T) - z(0) \sim \mathcal{N}(0, \sigma^2 T)$$

**In plain language:** The longer the branch (more evolutionary time), the more the trait is expected to change, and the variance of that change grows linearly with time. This is the defining property of Brownian motion.

The key insight of Felsenstein (1985) is that closely related species share more of their evolutionary history (more overlapping branches), making their trait values statistically correlated. The phylogenetic tree defines a covariance matrix among species, and this must be accounted for in any comparative analysis.

Extensions beyond BM include the Ornstein-Uhlenbeck process for modeling stabilizing selection (Hansen 1997; Butler & King 2004), Levy processes for modeling pulsed/jump evolution (Landis et al. 2013; Bastide & Didier 2023), and general stochastic diffusion models (Landis & Schraiber 2020). A comprehensive textbook treatment is available in Harmon (2019).

### 2.2 From Flat to Curved: Brownian Motion on Riemannian Manifolds

All of the classical phylogenetic BM models assume that traits evolve in flat Euclidean space. But what happens when the trait space is curved?

**In plain language:** Imagine ants doing a random walk on a flat table versus on the surface of a sphere. On the table, the random walk spreads out in all directions equally. On the sphere, the curvature constrains the walk — it can't go "through" the sphere, and distances along the surface differ from straight-line distances through the interior. The mathematics of random walks on curved surfaces is **Brownian motion on Riemannian manifolds**.

**Mathematically:** A Riemannian manifold *(M, g)* is a smooth space equipped with a **metric tensor** *g* that defines local notions of distance, angle, and volume. Brownian motion on *M* is the diffusion process whose infinitesimal generator is *(1/2) Delta_M*, where *Delta_M* is the **Laplace-Beltrami operator** — the generalization of the ordinary Laplacian to curved spaces.

In local coordinates *(x^1, ..., x^d)* with metric tensor *g_{ij}*, the Laplace-Beltrami operator is:

$$\Delta_M f = \frac{1}{\sqrt{\det g}} \, \partial_i \left( \sqrt{\det g} \, g^{ij} \, \partial_j f \right)$$

**In plain language:** This formula tells us how to compute "the average of a function over a tiny ball" on the curved surface. The metric tensor *g* appears because on a curved surface, "tiny balls" are not perfect circles — they're distorted by the curvature.

The transition density of this Brownian motion (the probability of moving from point *x* to point *y* in time *t*) is given by the **heat kernel** *p_t(x, y)*, which satisfies the heat equation on the manifold. For short times, the heat kernel has the asymptotic expansion (Minakshisundaram & Pleijel 1949):

$$p_t(x, y) \sim (4\pi t)^{-d/2} \exp\!\left(-\frac{d(x,y)^2}{4t}\right) \sum_{k=0}^{\infty} u_k(x,y) \, t^k$$

where *d(x, y)* is the geodesic distance and the coefficients *u_k* encode geometric information. The leading coefficient *u_0(x, x) = 1*, and the next coefficient *u_1(x, x) = R(x)/6* involves the **scalar curvature** *R(x)*, which measures how the volume of small balls deviates from flat space.

**In plain language:** The heat kernel says that for short times, Brownian motion on a manifold "looks like" Brownian motion in flat space (the Gaussian factor *exp(-d^2/4t)*), but with corrections that depend on the curvature. In regions of positive curvature (like a sphere), balls are smaller than in flat space, and the Brownian motion is slightly more concentrated.

The standard mathematical reference for stochastic analysis on manifolds is Hsu (2002).

### 2.3 Recent Work: Non-Euclidean Morphospaces

The recognition that morphological trait spaces are often non-Euclidean is gaining traction in evolutionary biology. Boyko & Rabosky (2025/2026) directly address this in their paper "The Geometry of Macroevolution: Phenotypic Evolution on Non-Euclidean Manifolds," developing methods for modeling phenotypic evolution within curved geometric frameworks. Boyko & Beaulieu (2025) use Lie group theory and Riemannian geometry to model evolution of the quantitative genetic covariance matrix (G-matrix) on phylogenies, where the space of positive-definite matrices is endowed with the Log-Euclidean metric.

Our approach is complementary: rather than working with a theoretically motivated manifold (like the space of covariance matrices), we use a **learned** manifold — the Riemannian structure induced by a VAE trained on empirical 3D beak shape data.

---

## 3. The VAE's Riemannian Geometry

### 3.1 How a VAE Decoder Induces a Riemannian Metric

The fundamental insight, developed by Arvanitidis, Hansen, & Hauberg (2018) in their paper "Latent Space Oddity," is that a VAE's decoder defines a natural Riemannian metric on the latent space.

**In plain language:** Imagine you have a rubber sheet (the latent space) that can be stretched and folded to cover a complex surface (the data space, i.e., the space of all possible beak shapes). The decoder is the mapping from the sheet to the surface. Where the decoder stretches the sheet a lot (small latent changes produce large data changes), the metric is large — meaning you "feel" more distance per unit of latent movement. Where the decoder compresses the sheet (latent changes produce small data changes), the metric is small.

**Mathematically:** Given a decoder function *f: R^d -> R^D* mapping latent codes *z* to data space, the **pullback metric** is:

$$M(z) = J(z)^T J(z)$$

where *J(z) = partial f / partial z* is the *D x d* Jacobian of the decoder. This *d x d* positive semi-definite matrix defines a Riemannian metric on the latent space.

**In plain language:** Each entry *M(z)_{ij}* tells you how much a small step in latent direction *i* combined with a step in direction *j* translates into distance in the decoded (beak shape) space. The Jacobian *J* captures how sensitive the decoder output is to changes in each latent dimension.

For stochastic decoders (as in VAEs with learned variance), the expected metric also includes a term from the variance function (Arvanitidis et al. 2018; Chadebec & Allassonniere 2022).

### 3.2 Our Data-Dependent Diagonal Metric

Computing the full pullback metric *J^T J* is expensive (it requires the full Jacobian of the decoder at every point). Instead, we use a **data-dependent diagonal metric** that approximates the same intuition using the VAE's posterior distribution.

**The formula** (implemented in `get_metric_tensor` in `R/.VAE_evo_model_v3_bayesian.R`):

$$G(z)_{jj} = \frac{1}{\displaystyle\sum_{i=1}^{N} \frac{1}{\sigma_{ij}^2} \exp\!\left(-\frac{d_{\text{Mah}}(z, \mu_i)^2}{\rho^2}\right) + \lambda}$$

where:

- *z* is a point in the 16-dimensional active latent space
- *mu_i* are the VAE posterior means for each of the *N* observed species (the "centroids")
- *sigma_{ij}^2* are the VAE posterior variances for species *i* in dimension *j*
- *d_Mah(z, mu_i)* is a Mahalanobis-like distance from *z* to centroid *mu_i*
- *rho* is a bandwidth parameter controlling the range of influence of each centroid
- *lambda* is a small regularization constant preventing division by zero

**Let's unpack this step by step:**

**Step 1: The Gaussian kernel weights.** The term *exp(-d_Mah(z, mu_i)^2 / rho^2)* is a Gaussian kernel centered at each species' centroid *mu_i*. It's large (close to 1) when *z* is near centroid *mu_i*, and decays to zero as *z* moves away. The parameter *rho* controls how quickly this decay happens.

**In plain language:** Each observed species "votes" on the local metric, and species that are closer to the point *z* have a stronger vote.

**Step 2: Weighting by inverse variance.** Each species' vote is weighted by *1/sigma_{ij}^2* — the inverse of its posterior variance in dimension *j*. Species with small variance (high certainty) in dimension *j* have more influence on the metric in that dimension.

**In plain language:** If the VAE is very confident about where a species sits in a particular latent dimension, that species has more say in defining the local geometry in that dimension.

**Step 3: The sum.** The denominator sums these weighted votes across all species. Where many species cluster closely (high data density), the sum is large, making *G(z)_{jj}* small.

**Step 4: The reciprocal.** Taking *1/(sum + lambda)* means:
- **In data-dense regions** (many nearby centroids, large sum): *G(z)_{jj}* is **small** — movement is "cheap"
- **In data-sparse regions** (few nearby centroids, small sum): *G(z)_{jj}* is **large** — movement is "expensive"

**In plain language:** The metric tensor acts like a "cost map" for movement through the latent space. Walking through regions full of real beak shapes is easy (low cost). Walking through empty regions — where no real beaks exist — is hard (high cost). This is exactly what we want for an evolutionary prior: ancestors should have had plausible beak shapes, not alien morphologies that no bird has ever had.

---

## 4. Path Energy vs Path Length

Two fundamental ways to measure the "size" of a curve on a Riemannian manifold play a central role in our formulation.

### 4.1 Definitions

Given a smooth curve *gamma: [0, 1] -> M* on a Riemannian manifold *(M, g)*, we define:

**Path length** (the total distance traveled along the curve):

$$L[\gamma] = \int_0^1 \sqrt{g_{\gamma(t)}\!\left(\gamma'(t), \gamma'(t)\right)} \, dt = \int_0^1 \|\gamma'(t)\|_g \, dt$$

**In plain language:** Imagine walking along the curve with a pedometer. The path length is the total distance your pedometer records, accounting for the manifold's curvature.

**Path energy** (the total "kinetic energy" along the curve):

$$E[\gamma] = \frac{1}{2} \int_0^1 g_{\gamma(t)}\!\left(\gamma'(t), \gamma'(t)\right) \, dt = \frac{1}{2} \int_0^1 \|\gamma'(t)\|_g^2 \, dt$$

**In plain language:** Instead of measuring total distance, the energy measures total "effort," with a quadratic penalty — going twice as fast costs four times as much energy.

### 4.2 The Cauchy-Schwarz Relationship

These two quantities are related by the **Cauchy-Schwarz inequality**:

$$L[\gamma]^2 = \left(\int_0^1 \|\gamma'(t)\|_g \cdot 1 \, dt\right)^2 \leq \left(\int_0^1 \|\gamma'(t)\|_g^2 \, dt\right) \cdot \left(\int_0^1 1^2 \, dt\right) = 2E[\gamma]$$

**In plain language:** The Cauchy-Schwarz inequality says that the square of the average is at most the average of the squares. Here, it means the squared path length is at most twice the path energy.

Equality holds if and only if *||gamma'(t)||_g* is constant — i.e., the curve has **constant speed**. For a constant-speed curve:

$$L[\gamma]^2 = 2E[\gamma]$$

### 4.3 Why Energy Is Preferred

Geodesics (the "straightest possible paths" on a manifold) are critical points of **both** functionals. But the energy functional is mathematically superior for several reasons:

1. **Unique parameterization:** Minimizing energy automatically produces constant-speed curves. Minimizing length doesn't — any reparameterization of a length-minimizer is also a length-minimizer. This ambiguity makes length harder to work with computationally.

2. **Nicer Euler-Lagrange equations:** The energy's Euler-Lagrange equation is the standard geodesic equation:

$$\ddot{\gamma}^k + \Gamma^k_{ij}(\gamma) \, \dot{\gamma}^i \, \dot{\gamma}^j = 0$$

where *Gamma^k_{ij}* are the Christoffel symbols. The length functional's Euler-Lagrange equation is more complicated and degenerate.

3. **Connection to Brownian motion:** As we'll see in the next section, the energy functional appears directly in the Onsager-Machlup action for Brownian motion. The length functional does not.

4. **Smoother optimization landscape:** The energy *||v||^2* is smooth everywhere (differentiable at *v = 0*), while the length *||v||* has a kink at *v = 0*, which can cause gradient issues in optimization.

**In plain language:** Using energy instead of length is like the difference between minimizing the sum of squared residuals (least squares) versus the sum of absolute residuals (least absolute deviations). Both find the same "best fit line," but least squares is smoother, has a unique solution, and connects to Gaussian probability theory.

---

## 5. The Onsager-Machlup Functional

### 5.1 Historical Context

The Onsager-Machlup (OM) functional was introduced by Lars Onsager and Stefan Machlup in their two-part 1953 paper on fluctuations in thermodynamic systems near equilibrium (Onsager & Machlup 1953; Machlup & Onsager 1953). Their key insight was that among all possible trajectories a stochastic system might follow, some are more probable than others, and the relative probability can be expressed through a variational principle — much like how classical mechanics describes the most likely trajectory of a particle through the principle of least action.

### 5.2 The Tube Probability Definition

**In plain language:** Imagine watching a pollen grain undergo Brownian motion in water. You draw a thin tube (of radius *epsilon*) around some reference curve in space. The question is: what is the probability that the pollen grain's actual trajectory stays inside this tube? The OM functional tells you how this probability depends on the shape of the reference curve.

**Mathematically:** Given a stochastic process *X_t* and two smooth reference curves *phi_1, phi_2*, the OM functional *L* is defined by:

$$\lim_{\epsilon \to 0} \frac{P\!\left(d(X_t, \phi_1(t)) \leq \epsilon \;\;\forall\, t \in [0,T]\right)}{P\!\left(d(X_t, \phi_2(t)) \leq \epsilon \;\;\forall\, t \in [0,T]\right)} = \exp\!\left(-\int_0^T L(\phi_1, \dot{\phi}_1) \, dt + \int_0^T L(\phi_2, \dot{\phi}_2) \, dt\right)$$

**In plain language:** As the tube gets thinner and thinner, the *ratio* of probabilities for two different reference curves converges to something determined by the OM functional *L*. A curve with a smaller OM action is more probable.

### 5.3 The OM Action in Flat Space

For a diffusion process in *R^d* satisfying the SDE

$$dX_t = b(X_t) \, dt + \sigma \, dW_t$$

with drift *b* and constant diffusion coefficient *sigma*, the OM Lagrangian is (Durr & Bach 1978; Takahashi & Watanabe 1981):

$$L(x, v) = \frac{1}{2\sigma^2} \|v - b(x)\|^2 + \frac{1}{2} \nabla \cdot b(x)$$

For a standard Wiener process (*b = 0*, *sigma = 1*), this simplifies to:

$$L(x, v) = \frac{1}{2} \|v\|^2$$

and the OM action is simply the **path energy**:

$$S_{\text{OM}}[\phi] = \frac{1}{2} \int_0^T \|\dot{\phi}(t)\|^2 \, dt = E[\phi]$$

**In plain language:** For ordinary Brownian motion in flat space (no drift, unit diffusion), the most probable path is the one that minimizes the total kinetic energy — which is a straight line traversed at constant speed (a geodesic in flat space). This makes intuitive sense: Brownian motion is "lazy" and prefers paths that require the least total effort.

### 5.4 The OM Action on Riemannian Manifolds

Now we come to the key result. On a Riemannian manifold *(M, g)*, Brownian motion is the diffusion with generator *(1/2) Delta_M*. The OM functional for this process was derived by Fujita & Kotani (1982):

**Theorem (Fujita & Kotani 1982).** The Onsager-Machlup Lagrangian for Brownian motion on a *d*-dimensional Riemannian manifold *(M, g)* is:

$$L(x, v) = \frac{1}{2} \|v\|_x^2 - \frac{1}{12} R(x)$$

where *||v||_x^2 = g_{ij}(x) v^i v^j* is the squared Riemannian norm of the velocity, and *R(x)* is the **scalar curvature** at *x*.

The OM action for a path *phi* is therefore:

$$S_{\text{OM}}[\phi] = \int_0^T \left[\frac{1}{2} \|\dot{\phi}(t)\|_{\phi(t)}^2 - \frac{1}{12} R(\phi(t))\right] dt$$

**Let's unpack each piece:**

**The energy term** *(1/2) ||v||^2*: This is identical to the flat-space case — it penalizes fast movement. The only difference is that the norm is now the *Riemannian* norm, so the cost of moving depends on the local metric.

**The scalar curvature correction** *-(1/12) R(x)*: This is the new ingredient. The scalar curvature *R(x)* measures how the volume of small geodesic balls deviates from flat space:

- *R(x) > 0* (positive curvature, like a sphere): small balls have *less* volume than in flat space. The *-R/12* term *reduces* the action, meaning paths through positively curved regions are *more* probable.
- *R(x) < 0* (negative curvature, like a saddle): small balls have *more* volume. The correction *increases* the action, making such paths *less* probable.

**In plain language:** Brownian motion on a curved surface is slightly biased by the curvature. In positively curved regions, the random walker is "squeezed" into a smaller volume and thus more concentrated — making paths through such regions more probable. The -(1/12)R term quantifies this bias exactly.

**Where does the -(1/12)R come from?** It arises from the short-time heat kernel expansion. The diagonal heat kernel has the expansion *p_t(x, x) ~ (4 pi t)^{-d/2} (1 + R(x) t/6 + ...)* (Minakshisundaram & Pleijel 1949; Vassilevich 2003). When this expansion is used to compute tube probabilities (which involve a conditioned or "bridge" process), the *R/6* coefficient transforms into *R/12* (see Hara & Takahashi 2016 for the full derivation).

### 5.5 The Most Probable Path

The path *phi** that maximizes the tube probability is the one that **minimizes** the OM action (Durr & Bach 1978):

$$\phi^* = \underset{\phi}{\text{argmin}} \; S_{\text{OM}}[\phi] = \underset{\phi}{\text{argmin}} \int_0^T \left[\frac{1}{2} \|\dot{\phi}\|^2 - \frac{1}{12} R(\phi)\right] dt$$

On a flat manifold (*R = 0*), this reduces to minimizing the path energy, and the solution is a geodesic. On a curved manifold, the curvature correction acts as a potential energy term that biases the most probable path toward regions of positive curvature.

**In plain language:** The most probable Brownian path is like a ball rolling on a landscape, where the "height" at each point is determined by *-R/12*. The ball prefers to roll through valleys (positive curvature) rather than over ridges (negative curvature), all while minimizing its total kinetic energy.

---

## 6. MAP Estimation: From Loss Function to Posterior

### 6.1 Finite-Dimensional MAP Estimation

In standard Bayesian inference, given data *D*, parameters *theta*, likelihood *p(D | theta)*, and prior *p(theta)*, the **posterior** is:

$$p(\theta \mid D) \propto p(D \mid \theta) \, p(\theta)$$

The **Maximum A Posteriori (MAP)** estimate is:

$$\theta_{\text{MAP}} = \underset{\theta}{\text{argmax}} \; p(\theta \mid D) = \underset{\theta}{\text{argmin}} \left[-\log p(D \mid \theta) - \log p(\theta)\right]$$

**In plain language:** The MAP estimate is the single most probable value of the parameters given the data. Finding it is equivalent to minimizing a sum of two terms: a "data fit" term (negative log-likelihood) and a "regularization" term (negative log-prior).

**Example:** If we have Gaussian data *Y ~ N(theta, sigma^2)* with a Gaussian prior *theta ~ N(0, tau^2)*, then:

$$-\log p(Y \mid \theta) = \frac{1}{2\sigma^2} \|Y - \theta\|^2 + \text{const}$$

$$-\log p(\theta) = \frac{1}{2\tau^2} \|\theta\|^2 + \text{const}$$

The MAP estimate minimizes their sum: *||Y - theta||^2 / (2 sigma^2) + ||theta||^2 / (2 tau^2)*. This is exactly **ridge regression** (L2-penalized least squares), with the penalty parameter *lambda = sigma^2 / tau^2* being the ratio of the observation noise variance to the prior variance.

### 6.2 Infinite-Dimensional MAP via the Onsager-Machlup Functional

In our problem, the "parameters" are not just numbers — they include entire *paths* through the latent space (one path per phylogenetic branch). This is an infinite-dimensional optimization problem, and the ordinary notion of a probability density doesn't apply (there's no Lebesgue measure in function spaces).

The rigorous solution, developed by Dashti et al. (2013) and further refined by Ayanbayev et al. (2022) and Kretschmann (2023), is that **the MAP estimate in function space is the minimizer of the Onsager-Machlup functional of the posterior measure**.

**In plain language:** In finite dimensions, the MAP is where the posterior density is highest. In infinite dimensions, there's no "density," but we can still ask: around which function do small balls have the highest posterior probability? The answer is the minimizer of the OM functional.

**Theorem (Dashti et al. 2013).** Let *mu* be a posterior measure on a separable Hilbert space, with a Gaussian prior whose Cameron-Martin space is *H_0*. Under mild conditions, the MAP estimator is:

$$\phi_{\text{MAP}} = \underset{\phi \in H_0}{\text{argmin}} \left[\frac{1}{2} \|\phi\|_{H_0}^2 + \Phi(\phi)\right]$$

where *||phi||_{H_0}* is the Cameron-Martin norm (encoding the prior) and *Phi* is the negative log-likelihood.

**In plain language:** This formula has the exact same structure as finite-dimensional MAP — it's still "minimize negative log-likelihood plus negative log-prior" — but the "prior penalty" is now the Cameron-Martin norm, which for Brownian motion is precisely the path energy.

### 6.3 The Full Bayesian Model for v3

We now show that the v3 loss function is the negative log-posterior of a well-defined Bayesian model. The loss function (from `R/.VAE_evo_model_v3_bayesian.R`) is:

$$\mathcal{L} = w_{\text{mani}} \cdot L_{\text{manifold}} + w_{\text{code}} \cdot L_{\text{code}} + w_{\text{troph}} \cdot L_{\text{trophic}} + w_{\text{tip}} \cdot L_{\text{tip}} + w_{\text{root}} \cdot L_{\text{root}}$$

We now derive the Bayesian interpretation of each term.

---

#### Term 1: Tip Loss (Gaussian Likelihood)

The tip loss measures how well the model's predicted tip latent codes match the observed ones:

$$L_{\text{tip}} = \frac{1}{N_{\text{tips}} \cdot d} \sum_{i=1}^{N_{\text{tips}}} \|Y_i - \hat{z}_i\|^2$$

where *Y_i* are the observed latent codes at the tips and *z_hat_i* are the predicted codes (from summing rates along the root-to-tip path).

**Bayesian interpretation:** This is the negative log-likelihood of a Gaussian observation model:

$$Y_i \mid z_i \sim \mathcal{N}(z_i, \sigma_{\text{tip}}^2 \, I)$$

Taking the negative log-likelihood:

$$-\log p(Y \mid z) = \frac{1}{2\sigma_{\text{tip}}^2} \sum_i \|Y_i - z_i\|^2 + \text{const}$$

Comparing with *w_tip * L_tip*, where *L_tip* is the MSE (sum of squares divided by the number of elements), we get:

$$w_{\text{tip}} = \frac{N_{\text{tips}} \cdot d}{2\sigma_{\text{tip}}^2} \quad \Longrightarrow \quad \sigma_{\text{tip}}^2 = \frac{N_{\text{tips}} \cdot d}{2 \, w_{\text{tip}}}$$

With *w_tip = 10*, *N_tips ≈ 2020*, and *d = 16*, this gives *sigma_tip^2 ≈ 1616*. This is a relatively **weak** likelihood (large observation variance), reflecting the fact that the latent codes have inherent uncertainty from the VAE encoding process.

**In plain language:** The tip loss says "the predicted tip values should be close to the observed ones, but we allow some slack." The weight *w_tip = 10* determines how much slack: higher weight means less slack (more precise observations).

---

#### Term 2: Manifold Energy (Onsager-Machlup Prior on Evolutionary Paths)

This is the central term. In v3, for each edge *e* with branch length *T_e*, we compute:

$$L_{\text{manifold}} = \text{mean}_{\text{edges, segs}} \left[\frac{\|\Delta z_k\|_{G(z_k)}^2}{T_e}\right]$$

where *Delta z_k* is the latent displacement over segment *k*, *||.||_G^2* is the squared Riemannian norm under the metric tensor *G*, and the mean is over all segments and all edges in the batch.

**Bayesian interpretation:** This corresponds to the Onsager-Machlup action for Brownian motion on the Riemannian manifold (neglecting the scalar curvature correction — see Section 7.4).

For a single edge *e*, the energy is:

$$E_e = \frac{1}{T_e} \int_0^1 \|\dot{z}_e(s)\|_{G(z_e(s))}^2 \, ds$$

where the integral is over the normalized path parameter *s in [0, 1]* and the *1/T_e* factor converts from normalized to real time (see derivation below).

**Derivation of the 1/T_e factor:** Let *s = t / T_e* be the normalized time parameter, so *t = T_e s* and *dt = T_e ds*. The real-time velocity is:

$$\frac{dz}{dt} = \frac{1}{T_e} \frac{dz}{ds}$$

The Onsager-Machlup energy in real time is:

$$E_e = \frac{1}{2} \int_0^{T_e} \left\|\frac{dz}{dt}\right\|_G^2 dt = \frac{1}{2} \int_0^1 \frac{1}{T_e^2} \left\|\frac{dz}{ds}\right\|_G^2 T_e \, ds = \frac{1}{2T_e} \int_0^1 \left\|\frac{dz}{ds}\right\|_G^2 ds$$

**In plain language:** The energy is the integral of squared velocity, divided by the branch length. This division is crucial for the Brownian motion interpretation: over longer branches, the same amount of total change is *less* surprising because there was more time for random fluctuations to accumulate. Mathematically, the variance of BM displacement grows linearly with time: *Var(z(T) - z(0)) = sigma^2 T*. So dividing the squared displacement by *T* normalizes for this expected scaling.

The discretized version (with *n_segs* grid points giving *n_segs - 1* segments of normalized step size *Delta s = 1/n_segs*) is:

$$E_e \approx \frac{1}{2T_e} \sum_{k=1}^{n_{\text{segs}}-1} \frac{\|\Delta z_k\|_G^2}{\Delta s} = \frac{n_{\text{segs}}}{2T_e} \sum_{k=1}^{n_{\text{segs}}-1} \|\Delta z_k\|_G^2$$

The negative log-prior for the evolutionary path on edge *e* is then:

$$-\log p(\text{path}_e) = \frac{E_e}{\sigma_{\text{evo}}^2} = \frac{1}{2\sigma_{\text{evo}}^2 T_e} \int_0^1 \left\|\frac{dz}{ds}\right\|_G^2 ds + \text{const}$$

The weight *w_mani* maps to the evolutionary rate precision: *w_mani ~ 1 / (2 sigma_evo^2)*. Larger *w_mani* means slower evolution (stronger constraint to stay near the ancestor), while smaller *w_mani* allows faster evolutionary change.

**In plain language:** This prior says "evolution is a random walk on the manifold, where the cost of moving through any region is determined by the local data density." In practice, this means ancestors are constrained to have had biologically plausible beak shapes (those supported by the VAE's training data), with the strength of this constraint controlled by the evolutionary rate parameter *sigma_evo^2*.

---

#### Term 3: Decoded Beak Shape Energy (Smoothness Prior in Phenotype Space)

$$L_{\text{code}} = \text{mean}_{\text{edges, segs}} \left[\frac{\|\Delta \text{dec}(z_k)\|^2}{T_e}\right]$$

where *dec(z_k)* is the VAE decoder applied to the latent path points, outputting the predicted beak shape codes.

**Bayesian interpretation:** This is an *additional* Brownian motion prior, but in the decoded phenotype space rather than the latent space:

$$-\log p_{\text{code}}(\text{path}_e) = \frac{\alpha_{\text{code}}}{T_e} \int_0^1 \left\|\frac{d}{ds}\text{dec}(z_e(s))\right\|^2 ds$$

**In plain language:** The manifold energy already constrains paths to stay in plausible regions of latent space. The decoded energy adds a complementary constraint: the *actual beak shapes* along the path should also change smoothly. This prevents the optimizer from finding paths that are smooth in latent space but produce jerky changes in the decoded beak morphology.

---

#### Term 4: Decoded Trophic Energy (Smoothness Prior in Ecological Niche Space)

$$L_{\text{trophic}} = \text{mean}_{\text{edges, segs}} \left[\frac{\|\Delta \text{softmax}(\text{troph}(z_k))\|^2}{T_e}\right]$$

where *troph(z_k)* is the trophic niche logit output of the decoder, and softmax converts these to probabilities.

**Bayesian interpretation:** A Brownian motion smoothness prior on the trophic niche probabilities:

$$-\log p_{\text{troph}}(\text{path}_e) = \frac{\alpha_{\text{troph}}}{T_e} \int_0^1 \left\|\frac{d}{ds}\text{softmax}(\text{troph}(z_e(s)))\right\|^2 ds$$

**In plain language:** Trophic niche (diet type) should change gradually along evolutionary branches, not jump abruptly. A lineage doesn't typically switch from "frugivore" to "piscivore" and back in a short time. This prior penalizes rapid switches in the predicted trophic probabilities.

---

#### Term 5: Root Loss (Gaussian Prior on Root State)

$$L_{\text{root}} = \|\mathbf{r}\|^2 = \sum_{j=1}^{d} r_j^2$$

where **r** is the learnable root state vector.

**Bayesian interpretation:** A Gaussian prior centered at the origin:

$$\mathbf{r} \sim \mathcal{N}(\mathbf{0}, \sigma_{\text{root}}^2 \, I)$$

$$-\log p(\mathbf{r}) = \frac{1}{2\sigma_{\text{root}}^2} \|\mathbf{r}\|^2 + \text{const}$$

With weight *w_root = 1/100*, we get *sigma_root^2 = 50*.

**In plain language:** The root of the tree (the common ancestor of all birds in the dataset) should have a latent code near the origin of the latent space. The prior is weak (*sigma_root^2 = 50* is quite large), reflecting our genuine uncertainty about the ancestral beak form.

---

### 6.4 The Complete Bayesian Model

Putting it all together, the complete posterior is:

$$p(\theta \mid Y) \propto \underbrace{p(Y \mid z_{\text{tips}})}_{\text{tip likelihood}} \cdot \underbrace{p(\mathbf{r})}_{\text{root prior}} \cdot \prod_{e \in \text{edges}} \underbrace{p_{\text{mani}}(\text{path}_e)}_{\text{manifold BM prior}} \cdot \underbrace{p_{\text{code}}(\text{path}_e)}_{\text{phenotype smoothness}} \cdot \underbrace{p_{\text{troph}}(\text{path}_e)}_{\text{trophic smoothness}}$$

where *theta = {rates, a, b, r}* are the learnable parameters (evolutionary rates, path curvature coefficients, and root state).

The negative log-posterior is:

$$-\log p(\theta \mid Y) = \frac{1}{2\sigma_{\text{tip}}^2} \sum_i \|Y_i - \hat{z}_i\|^2 + \frac{1}{2\sigma_{\text{root}}^2}\|\mathbf{r}\|^2 + \sum_e \left[\frac{E_e^{\text{mani}}}{\sigma_{\text{evo}}^2} + \frac{\alpha_{\text{code}} \, E_e^{\text{code}}}{1} + \frac{\alpha_{\text{troph}} \, E_e^{\text{troph}}}{1}\right] + \text{const}$$

Minimizing this is exactly what `optim_adam` does when it minimizes the v3 loss function.

### 6.5 Loss Weights as Precision Parameters

| Loss weight | Value | Bayesian parameter | Interpretation |
|---|---|---|---|
| `tip_weight` | 10 | *1 / (2 sigma_tip^2)* (up to constants) | Observation precision |
| `manifold_weight` | 1.0 | *1 / (2 sigma_evo^2)* | Evolutionary rate precision (inverse rate) |
| `code_weight` | 1/64 | *alpha_code* | Phenotype smoothness precision |
| `trophic_weight` | 1/10 | *alpha_trophic* | Trophic niche smoothness precision |
| `root_weight` | 1/100 | *1 / (2 sigma_root^2)* | Root state precision |

**In plain language:** Each weight is like a dial controlling the strength of a belief:
- High `tip_weight` = "I trust the observed data strongly"
- High `manifold_weight` = "I believe evolution is slow (small rate)"
- High `code_weight` = "I believe beak shape changes gradually"
- High `trophic_weight` = "I believe diet changes gradually"
- High `root_weight` = "I believe the ancestor was near the latent space origin"

---

## 7. What Changed from v2 to v3

### 7.1 Path Length to Path Energy

The single most important change is in `get_manifold_dist` / `get_manifold_energy`.

**v2** (`R/.VAE_evo_model_v2.0.R`, line 274):
```r
get_manifold_dist <- function(vel, metric) {
  (vel*metric*vel)$sum(dim = 2)$sqrt()
}
```

This computes the **Riemannian norm** *||vel||_G = sqrt(vel^T G vel)*, giving the **path length** integrand.

**v3** (`R/.VAE_evo_model_v3_bayesian.R`, line 270):
```r
get_manifold_energy <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)
}
```

This computes the **squared Riemannian norm** *||vel||_G^2 = vel^T G vel*, giving the **path energy** integrand.

**The probabilistic consequence:**

- **Path length** (v2): *exp(-L/tau)* is a **Laplace-type** prior. The penalty is linear in displacement magnitude, corresponding to heavy-tailed (Laplace) distributions. This allows occasional large jumps more readily than a Gaussian model.

- **Path energy** (v3): *exp(-E/sigma^2)* is a **Gaussian** prior. The penalty is quadratic in displacement, corresponding to the true transition density of Brownian motion. Large jumps are exponentially more penalized.

**In plain language:** v2 was like using absolute error (|x|); v3 uses squared error (x^2). Both penalize deviations, but squared error penalizes large deviations much more strongly. The squared version is the one that corresponds to actual Brownian motion.

### 7.2 Branch Length Scaling

**v2** divided distances by `blens / n_segs` (the real time per segment), computing instantaneous velocity:
```r
dists <- get_manifold_dist(zs[[3]], met) / (x[[5]]$unsqueeze(-1) / self$n_segs)
```

**v3** divides energies by `blens` (the full branch length), computing the energy with proper *1/T_e* BM scaling:
```r
manifold_energies <- get_manifold_energy(zs[[3]], met) / blens_expanded
```

The *1/T_e* factor is essential for the Brownian motion interpretation (Section 6.3): longer branches should allow more evolutionary change.

### 7.3 Explicit Loss Weights

v3 introduces named, documented loss weights that map directly to Bayesian parameters:
```r
manifold_weight <- 1.0
code_weight <- 1/64
trophic_weight <- 1/10
tip_weight <- 10
root_weight <- 1/100
```

### 7.4 The Scalar Curvature Correction: Why It's Omitted

The full Onsager-Machlup action includes a *-(1/12) R(z)* scalar curvature term (Section 5.4). We omit it for three reasons:

1. **Computational cost:** Computing the scalar curvature requires second derivatives of the metric tensor, which requires differentiating through the Gaussian kernel sums — expensive on GPU.

2. **Magnitude:** For our diagonal metric, the curvature corrections are expected to be small relative to the energy term, especially since our metric is a smooth kernel density estimate rather than a sharply curved surface.

3. **Standard practice:** The curvature correction is routinely omitted in practical applications of the OM functional (Dashti et al. 2013; Arvanitidis et al. 2018). In the statistics literature on infinite-dimensional MAP estimation, the curvature correction appears as a higher-order term that doesn't affect the MAP estimator to leading order.

This omission means our model corresponds to the OM action for BM on a manifold *as if the manifold were locally flat* — i.e., we use the leading-order approximation. This is increasingly accurate as the metric varies slowly relative to the step size.

---

## 8. The Rho Annealing Schedule as Posterior Tempering

The bandwidth parameter *rho* in the metric tensor controls how sharply the metric varies across the latent space. The training procedure anneals *rho* from a large initial value (*3 x max_min_dist*) to a small target value (*max_min_dist / 3*) using a cosine schedule.

**Bayesian interpretation:** This is a form of **posterior tempering** (also known as simulated annealing), a well-established technique in Bayesian computation.

Define a tempered posterior:

$$p_\beta(\theta \mid Y) \propto p(Y \mid \theta) \, p(\theta)^\beta$$

where *beta in [0, 1]* is the inverse temperature. At *beta = 0*, we have a flat (improper) prior — the posterior reduces to the likelihood. At *beta = 1*, we have the full posterior.

In our case, the annealing of *rho* modifies the metric tensor, which plays the role of the prior. When *rho* is large:

- The Gaussian kernels *exp(-d^2/rho^2)* overlap heavily
- All centroids contribute roughly equally everywhere
- The metric *G(z)* is approximately constant (Euclidean)
- The manifold prior is nearly **flat** (uninformative)

When *rho* is small:

- Only nearby centroids contribute
- The metric varies sharply — high in sparse regions, low in dense regions
- The manifold prior is strongly **informative**, forcing paths through data-dense regions

**In plain language:** We start training with a nearly flat (Euclidean) geometry and gradually "turn on" the manifold structure. This prevents the optimizer from getting stuck in bad local optima early on (when the loss landscape would be very rough with the full manifold), similar to how simulated annealing in metallurgy starts at high temperature (allowing broad exploration) and gradually cools to find the global optimum.

---

## 9. Gradient Accumulation for Optimization

### 9.1 Why Mini-Batching Is Needed

The metric tensor computation (Section 3.2) requires computing the Mahalanobis distance from every path midpoint to every centroid — an *O(batch * N_species * d * n_segs)* operation. With ~4,000 edges, ~2,000 species, 16 dimensions, and 50 segments, doing this for all edges simultaneously would exceed GPU memory. The solution is to process edges in mini-batches.

### 9.2 The Gradient Accumulation Strategy

The training loop in the v3 code (and v2 before it) uses **gradient accumulation**: gradients from all mini-batches are summed before a single optimizer step per epoch.

```r
for(epoch in 1:n_epoch) {
  optim1$zero_grad()            # (1) Zero gradients ONCE per epoch
  coro::loop(for (b in bill_dl) {
    res <- mod(b)
    loss <- ...                  # (2) Compute loss for this batch
    loss$backward()              # (3) ACCUMULATE gradients (adds to existing)
  })
  optim1$step()                  # (4) ONE optimizer step per epoch
}
```

**This is NOT mini-batch SGD.** In standard SGD, you update parameters after *each* batch. Here, parameters are held **fixed** across all batches within an epoch, and gradients are accumulated. This is full-batch gradient descent with a memory-saving trick.

### 9.3 Mathematical Proof of Equivalence

**Claim:** The accumulated gradient equals a constant times the full-batch gradient.

**Setup:** Let the total loss over all *N* edges decompose into *K* batches *B_1, ..., B_K* of size *n* (assuming equal batch sizes for simplicity). The loss for batch *B_k* is:

$$\mathcal{L}_k(\theta) = \frac{1}{|B_k|} \sum_{e \in B_k} \ell_e(\theta) + w_{\text{tip}} \cdot L_{\text{tip}}(\theta) + w_{\text{root}} \cdot L_{\text{root}}(\theta)$$

where *l_e* is the per-edge distance/energy loss, and *L_tip* and *L_root* are computed on the full dataset in every batch.

**Step 1: Gradient of the per-edge terms.** The accumulated gradient of the per-edge (dist/code/trophic) terms is:

$$\sum_{k=1}^{K} \nabla_\theta \left[\frac{1}{n} \sum_{e \in B_k} \ell_e(\theta)\right] = \frac{1}{n} \sum_{e=1}^{N} \nabla_\theta \, \ell_e(\theta) = \frac{N}{n} \cdot \frac{1}{N} \sum_{e=1}^{N} \nabla_\theta \, \ell_e(\theta) = K \cdot \nabla_\theta \, \mathcal{L}_{\text{full}}^{\text{edges}}$$

**In plain language:** Since every edge appears in exactly one batch and we're averaging within each batch, summing these averages across *K* batches gives *K* times the overall average. The key step is that the batches *partition* the full dataset (no overlap, no gaps), so summing partial averages gives a scaled version of the grand average.

**Step 2: Gradient of the tip and root terms.** These are computed identically in every batch (they don't depend on which edges are in the batch), so:

$$\sum_{k=1}^{K} \nabla_\theta \, L_{\text{tip}}(\theta) = K \cdot \nabla_\theta \, L_{\text{tip}}(\theta)$$

and similarly for the root loss.

**Step 3: The total accumulated gradient.** Combining:

$$\text{accumulated gradient} = K \cdot \nabla_\theta \, \mathcal{L}_{\text{full}}(\theta)$$

Since *K* is a constant, the accumulated gradient points in **exactly the same direction** as the full-batch gradient and differs only by a constant factor. This factor effectively scales the learning rate by *K* — which is fine, since the learning rate is a tunable hyperparameter.

**In plain language:** If you split a restaurant bill equally among friends, each person computes their share. Adding up all the shares gives the total bill. Similarly, each mini-batch computes its "share" of the gradient, and adding them up recovers the full gradient (times the number of batches).

### 9.4 Formal Statement

**Proposition.** Let *theta* be fixed model parameters. If the *N* edges are partitioned into *K* disjoint batches of equal size *N/K*, and gradients are computed for each batch without updating *theta* between batches, then:

$$\sum_{k=1}^{K} \nabla_\theta \mathcal{L}_k(\theta) = K \cdot \nabla_\theta \mathcal{L}_{\text{full}}(\theta)$$

This holds **exactly** (not in expectation) due to the linearity of the gradient operator and the disjoint partition of edges.

The equivalence relies on one critical assumption: **model parameters are not updated between batches**. If they were, the gradients would be evaluated at different points in parameter space, breaking the exact equivalence and reducing to standard mini-batch SGD. This is precisely why `optim1$zero_grad()` is called *once* at the beginning of the epoch (not after each batch), and `optim1$step()` is called *once* at the end.

**Key references:** The linearity property underlying gradient accumulation is a direct consequence of the chain rule and the additivity of derivatives. For formal treatments of stochastic optimization more broadly, see Bottou, Curtis, & Nocedal (2018) and Robbins & Monro (1951). For the specific technique of gradient accumulation, see Kingma & Ba (2015) where it is implicitly used in the description of Adam.

### 9.5 Interaction with the Adam Optimizer

Adam (Kingma & Ba 2015) maintains running averages of the first and second moments of the gradient:

$$m_t = \beta_1 \, m_{t-1} + (1 - \beta_1) \, g_t$$
$$v_t = \beta_2 \, v_{t-1} + (1 - \beta_2) \, g_t^2$$

with bias-corrected estimates *m_hat = m / (1 - beta_1^t)* and *v_hat = v / (1 - beta_2^t)*, and the update:

$$\theta_{t+1} = \theta_t - \eta \, \frac{\hat{m}_t}{\sqrt{\hat{v}_t} + \epsilon}$$

Since Adam normalizes by the second moment, the constant factor of *K* in the accumulated gradient is largely absorbed:

$$\frac{K \cdot g}{\sqrt{K^2 \cdot g^2 + \epsilon}} \approx \frac{g}{\sqrt{g^2 + \epsilon/K^2}} \approx \text{sign}(g)$$

for sufficiently large gradients. In practice, the *K* factor slightly affects the effective step size but not the direction, so the optimizer converges to the same minimum.

### 9.6 Observations on the Implementation

**Redundant tip loss computation:** The tip loss uses the full tip matrix and full tip data in every batch, not just the batch edges. This means the same tip loss gradient is computed *K* times per epoch. An optimization would be to compute it once outside the batch loop, but this doesn't affect correctness — it just wastes some GPU cycles.

**Last batch size:** With `drop_last = FALSE` and batch size 104 over ~4,039 edges, there are 38 full batches and one last batch of ~87 edges. Edges in the last batch are slightly upweighted (by factor 104/87 ≈ 1.20) in the per-edge loss terms. This is a standard and negligible artifact of unequal batch sizes.

**The a/b parameters:** The path curvature parameters `self$a` and `self$b` are indexed by batch (`self$a[x[[6]]]`). Since each edge appears in exactly one batch per epoch, the accumulated gradient correctly fills in all rows of `self$a$grad` and `self$b$grad` with no overlap or gaps.

---

## 10. Summary and Open Questions

### 10.1 What the Bayesian Interpretation Buys Us

1. **Interpretable priors.** Each loss term now has a clear probabilistic meaning. The manifold energy is a Brownian motion prior on evolutionary paths. The decoder losses are smoothness priors in phenotype and ecological niche space. The tip loss is a Gaussian likelihood. The root loss is a Gaussian prior.

2. **Principled weight selection.** The loss weights are no longer arbitrary tuning parameters — they are precision (inverse-variance) parameters with interpretable units. In principle, they could be estimated from the data (empirical Bayes) or given informative values based on domain knowledge about evolutionary rates.

3. **Path to full Bayesian inference.** With a proper posterior defined, one could in principle go beyond MAP estimation to full posterior sampling using techniques like Hamiltonian Monte Carlo (HMC) or its Riemannian generalization (Girolami & Calderhead 2011). This would provide uncertainty quantification for ancestral state reconstructions.

### 10.2 The Scalar Curvature Correction

The full OM action includes a *-(1/12) R(x)* curvature term that we currently neglect (Section 7.4). Computing this for our kernel-based metric would require:

1. Second derivatives of the metric tensor (Hessians of the Gaussian kernel sums)
2. Christoffel symbols and their derivatives
3. Contraction to the Riemann tensor and then the scalar curvature

This is computationally feasible (all operations are differentiable and could be done in torch) but expensive. Whether the correction materially changes the MAP estimate is an empirical question worth investigating.

### 10.3 Possible Extensions

- **Empirical Bayes for loss weights:** The precision parameters could be estimated by maximizing the marginal likelihood (integrating out the paths). This would provide data-driven weight selection.

- **Variable evolutionary rates:** Currently *sigma_evo^2* is a single global parameter. A natural extension would be to allow rates to vary across branches or clades (analogous to relaxed clock models in molecular evolution).

- **Alternative priors:** The Gaussian (BM) prior on paths could be replaced by Ornstein-Uhlenbeck (incorporating selection toward an optimum in latent space) or Levy processes (allowing punctuated jumps).

---

## 11. References

### Original Onsager-Machlup Papers

- Onsager, L. and Machlup, S. (1953). "Fluctuations and Irreversible Processes." *Physical Review*, 91(6): 1505-1512.
- Machlup, S. and Onsager, L. (1953). "Fluctuations and Irreversible Processes. II. Systems with Kinetic Energy." *Physical Review*, 91(6): 1512-1515.

### Onsager-Machlup Theory and Most Probable Paths

- Durr, D. and Bach, A. (1978). "The Onsager-Machlup function as Lagrangian for the most probable path of a diffusion process." *Communications in Mathematical Physics*, 60: 153-170.
- Takahashi, Y. and Watanabe, S. (1981). "The probability functionals (Onsager-Machlup functions) of diffusion processes." In *Stochastic Integrals*, Lecture Notes in Mathematics 851, Springer.
- Fujita, T. and Kotani, S. (1982). "The Onsager-Machlup function for diffusion processes." *Journal of Mathematics of Kyoto University*, 22: 115-130.
- Hara, K. and Takahashi, Y. (2016). "Stochastic analysis in a tubular neighborhood or Onsager-Machlup functions revisited." [arXiv:1610.06670](https://arxiv.org/abs/1610.06670).

### Heat Kernel Asymptotics

- Minakshisundaram, S. and Pleijel, A. (1949). "Some properties of the eigenfunctions of the Laplace operator on Riemannian manifolds." *Canadian Journal of Mathematics*, 1: 242-256.
- Vassilevich, D.V. (2003). "Heat kernel expansion: user's manual." *Physics Reports*, 388(5-6): 279-360. [arXiv:hep-th/0306138](https://arxiv.org/abs/hep-th/0306138).

### Stochastic Analysis on Manifolds

- Hsu, E.P. (2002). *Stochastic Analysis on Manifolds*. Graduate Studies in Mathematics 38, AMS.

### MAP Estimation in Infinite Dimensions

- Dashti, M., Law, K.J.H., Stuart, A.M., and Voss, J. (2013). "MAP estimators and their consistency in Bayesian nonparametric inverse problems." *Inverse Problems*, 29(9): 095017. [arXiv:1303.4795](https://arxiv.org/abs/1303.4795).
- Ayanbayev, B., Klebanov, I., Lie, H.C., and Sullivan, T.J. (2022). "Gamma-convergence of Onsager-Machlup functionals. Part I: With applications to maximum a posteriori estimation in Bayesian inverse problems." *Inverse Problems*. [arXiv:2108.04597](https://arxiv.org/abs/2108.04597).
- Kretschmann, R. (2023). "Are minimizers of the Onsager-Machlup functional strong posterior modes?" *SIAM/ASA Journal on Uncertainty Quantification*, 11(4). [arXiv:2212.04275](https://arxiv.org/abs/2212.04275).

### VAE Riemannian Geometry

- Arvanitidis, G., Hansen, L.K., and Hauberg, S. (2018). "Latent Space Oddity: on the Curvature of Deep Generative Models." *International Conference on Learning Representations (ICLR)*. [arXiv:1710.11379](https://arxiv.org/abs/1710.11379).
- Shao, H., Kumar, A., and Fletcher, P.T. (2018). "The Riemannian Geometry of Deep Generative Models." *CVPR Workshop*.
- Chen, N., Klushyn, A., Kurle, R., Jiang, X., Bayer, J., and van der Smagt, P. (2020). "Learning Flat Latent Manifolds with VAEs." *ICML 2020*.
- Chadebec, C. and Allassonniere, S. (2022). "A Geometric Perspective on Variational Autoencoders." *NeurIPS 2022*. [arXiv:2209.07370](https://arxiv.org/abs/2209.07370).

### Phylogenetic Comparative Methods

- Felsenstein, J. (1985). "Phylogenies and the Comparative Method." *The American Naturalist*, 125(1): 1-15.
- Hansen, T.F. (1997). "Stabilizing Selection and the Comparative Analysis of Adaptation." *Evolution*, 51(5): 1341-1351.
- Butler, M.A. and King, A.A. (2004). "Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution." *The American Naturalist*, 164(6): 683-695.
- Pagel, M. (1999). "Inferring the Historical Patterns of Biological Evolution." *Nature*, 401(6756): 877-884.
- Landis, M.J., Schraiber, J.G., and Liang, M. (2013). "Phylogenetic Analysis Using Levy Processes." *Systematic Biology*, 62(2): 193-204.
- Bastide, P. and Didier, G. (2023). "The Cauchy Process on Phylogenies." *Systematic Biology*, 72(6): 1296-1315.
- Landis, M.J. and Schraiber, J.G. (2020). "Beyond Brownian Motion and the Ornstein-Uhlenbeck Process." *The American Naturalist*, 195(2): 145-165.
- Harmon, L.J. (2019). *Phylogenetic Comparative Methods: Learning from Trees*. [Online textbook](https://lukejharmon.github.io/pcm/).

### Non-Euclidean Evolutionary Models

- Boyko, J.D. and Rabosky, D.L. (2025/2026). "The Geometry of Macroevolution: Phenotypic Evolution on Non-Euclidean Manifolds." *The American Naturalist*. [DOI:10.1086/740145](https://doi.org/10.1086/740145).
- Boyko, J.D. and Beaulieu, J.M. (2025). "Multivariate Trait Evolution: Models for the Evolution of the Quantitative Genetic G-matrix on Phylogenies." *Evolution Letters*, 9(6): 706-716.

### The Published DeepBills Work

- Dinnage, R. and Kleineberg, M. (2025). "Generative AI Extracts Ecological Meaning from the Complex Three Dimensional Shapes of Bird Bills." *PLoS Computational Biology*, 21(3): e1012887. [DOI:10.1371/journal.pcbi.1012887](https://doi.org/10.1371/journal.pcbi.1012887).
- Park, J.J., Florence, P., Straub, J., Newcombe, R., and Lovegrove, S. (2019). "DeepSDF: Learning Continuous Signed Distance Functions for Shape Representation." *CVPR*, 165-174.

### Riemannian Geometry Textbooks

- do Carmo, M.P. (1992). *Riemannian Geometry*. Birkhauser.
- Gallot, S., Hulin, D., and Lafontaine, J. (2004). *Riemannian Geometry*, 3rd ed. Springer.

### Optimization

- Kingma, D.P. and Ba, J. (2015). "Adam: A Method for Stochastic Optimization." *ICLR*. [arXiv:1412.6980](https://arxiv.org/abs/1412.6980).
- Bottou, L., Curtis, F.E., and Nocedal, J. (2018). "Optimization Methods for Large-Scale Machine Learning." *SIAM Review*, 60(2): 223-311. [arXiv:1606.04838](https://arxiv.org/abs/1606.04838).
- Robbins, H. and Monro, S. (1951). "A Stochastic Approximation Method." *Annals of Mathematical Statistics*, 22(3): 400-407.
- Girolami, M. and Calderhead, B. (2011). "Riemann manifold Langevin and Hamiltonian Monte Carlo methods." *Journal of the Royal Statistical Society: Series B*, 73(2): 123-214.

### Evolution in Latent Spaces

- Ding, D. et al. (2019). "Deciphering Protein Evolution and Fitness Landscapes with Latent Space Models." *Nature Communications*, 10: 5644.
- Xie, T., Richman, H., Gao, J., Matsen, F.A. IV, and Zhang, C. (2025). "PhyloVAE: Unsupervised Learning of Phylogenetic Trees via Variational Autoencoders." *ICLR 2025*. [arXiv:2502.04730](https://arxiv.org/abs/2502.04730).
