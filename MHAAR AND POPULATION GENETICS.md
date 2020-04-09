# MHAAR AND POPULATION GENETICS

$$
\newcommand{\v}[1]{\boldsymbol{\mathbf{#1}}}
\newcommand{\nc}[2]{\newcommand{#1}{#2}}
\nc{\vv}{\v{v}}
$$

### MHAAR

Let $\eta(\theta)$ be a **prior** and let $z_1, \ldots, z_T$ be **latent variables**, denoted collectively as $z=z_{1:T}$.
$$
\pi(\theta, z) \propto \eta(\theta)\prod_{t=1}^T \gamma_{t, \theta}(z_t)
$$
For every auxiliary variable $z_t$ sample $M-1$ **auxiliary variables** candidates $u_t^{(1)}, \ldots, u_t^{(M-1)}$ from a **proposal distribution** $q_{t, \theta, \theta'}$. We define
$$
\begin{align}
\vv = v_{1:T}^{(1:M)} &= (v_1^{(1)}, v_{1}^{(2)}, \ldots, v_{1}^{(M)}, \ldots \ldots, v_T^{(1)}, v_T^{(2)}, \ldots, v_T^{(M)}) \\
&=(z_1, u_1^{(1)}, \ldots, u_1^{(M-1)}, \ldots\ldots, z_T, u_T^{(2)}, \ldots,u_T^{(M-1)})
\end{align}
$$
Define **weights**
$$
w_{t, \theta, \theta'}(v) = \frac{\gamma_{t,\theta}(v)}{q_{t,\theta,\theta'}(v)}
$$
when there is **no annealing** the acceptance ratio becomes
$$
\mathring{r}_{\vv}(\theta, \theta') = \frac{q(\theta', \theta)\eta(\theta')}{q(\theta, \theta')\eta(\theta)} \prod_{t=1}^T \frac{\sum_{i=1}^M\frac{\gamma_{t,\theta'}(v^{(i)})}{q_{t, \theta, \theta'}(v^{(i)})}}{\sum_{i=1}^M \frac{\gamma_{t, \theta}(v^{(j)})}{q_{t,\theta, \theta'}(v^{(j)})}}
$$

### POPULATION GENETICS

In our case the parameter of interest is $\phi$ and the likelihood for a single nucleotide position $x_j$ is given by
$$
p(x_j\mid \phi) = \int p(x_j, H_j\mid \phi)d H_j = \int p(x_j\mid H_j, \phi) p(H_j\mid \phi) dH_j = \int p(H_j\mid\phi) dH_j
$$
assuming independence the full likelihood can be written as
$$
p(\v{x}\mid \phi) = \prod_{i=1}^Tp(x_j\mid \phi) = \prod_{t=1}^T\int p(H_j\mid \phi) dH_j
$$
Suppose we have a bunch of nucleotide positions