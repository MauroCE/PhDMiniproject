$$
\newcommand{\v}[1]{\boldsymbol{\mathbf{#1}}}
\newcommand{\nc}[2]{\newcommand{#1}{#2}}
\nc{\vv}{\v{v}}
\nc{\normal}{\mathcal{N}}
\newcommand{\vphi}{\v{\phi}}
\newcommand{\cH}{\mathcal{H}}
$$

# Population Genetics

Parameter of interest
$$
\vphi = (\log\theta, \log r, \log t_f)^\top
$$
Assume a multivariate normal prior on $\vphi$
$$
p(\vphi) =\normal(\v{0}, I_3; \vphi)
$$
Our data is composed of $100$ files
$$
D = \left\{D_1, \ldots, D_{100}\right\} = \left\{\texttt{gtree_file1}, \ldots, \texttt{gtree_file100}\right\}
$$
where each file contains the sequence types and their multiplicity (frequency) for a single genetic locus on a specific chromosome. For instance, this is $\texttt{gtree_file1}$ for a stable population

```sh
0 4 : 000001000000011000000
0 1 : 000001000000111000000
0 1 : 000001000100000000100
0 2 : 000001001000011000000
0 1 : 000001010011001100000
0 4 : 000101110010001000000
0 1 : 000101110010001000010
0 2 : 001001000000001010000
0 1 : 010011010010001000000
0 3 : 100000000000000001001
```

This represents the current state of the sample, at a specific locus. `0` means that the nucleotide is ancestral, while `1` means that there has been a mutation. We assume that there cannot be back-mutations (i.e. `0` mutates, becoming `1` and then mutates back to `0`). 

Our aim is to compute the posterior
$$
p(\vphi\mid D) \propto p(\vphi) p(D \mid \vphi) = p(\vphi) \prod_{d=1}^{100} p(D_d\mid \vphi)
$$
Each likelihood term $p(D_d\mid \vphi)$ is the marginal distribution of a joint distribution
$$
p(D_d\mid \vphi) = \int_{\cH_d} p(D_d, H_d \mid \vphi)dH_d = \int_{\mathcal{H}_d} p(D_d \mid H_d, \vphi) p(H_d \mid \vphi) d H_d = \int_{\cH_d} p(H_d\mid \vphi) dH_d
$$
where $H_d$ is a possible history of coalescence 
$$
H_d :=\left\{H_d^{(0)}, H_{d}^{(-1)}, \ldots, H_{d}^{(-m)}\right\}
$$
and each $H_d^{(j)}$ is an ancestral configuration at some event (mutation, or coalescence). So that $H_d^{(0)}$ is the current state, i.e. the state of the sample, and $H_{d}^{(-m)}$ is the configuration of the MRCA. By ancestral configuration we mean that $H_{d}^{(j)}$ is an unordered list of sequence types available at time $t=j$.
$$
p(\vphi\mid D) \propto p(\vphi) \prod_{d=1}^{100} p(D_d\mid \vphi)
$$

$$
p(D_d\mid \vphi) \approx \frac{1}{\text{isss}} \sum_{i=1}^{\text{isss}} \left[\frac{p(H_d^{(0), i}\mid H_d^{(-1),i})}{q(H_d^{(-1),i}\mid H_d^{(0),i})}\cdots \frac{p(H_d^{(-m+1),i}\mid H_d^{(-m),i})}{q(H_d^{(-m),i}\mid H_d^{(-m+1),i})}p(H_d^{(-m),i})\right]
$$

$$
H_d^{(j), k}
$$


$$
H_d^{(0)}
$$

$$
H_d^{(-1)}, \ldots, H_d^{(-m)}
$$

$$
\theta^\prime \theta'
$$



# New

Dataset
$$
D = \left\{\texttt{gtree_file1}, \ldots, \texttt{gtree_file100}\right\}
$$
Prior
$$
p(\phi) = \normal(\v{0}, I_3; \phi)
$$
Likelihood
$$

$$


