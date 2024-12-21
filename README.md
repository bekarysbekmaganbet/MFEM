## Solving the Allen-Cahn equation with MFEM

### Mathematical description
Let $\Omega$ be a bounded domain in $\mathbb{R}^2$ or $\mathbb{R}^3$. Let $\varepsilon > 0$ be a fixed small positive constant.

Consider the functional $J: H^1(\Omega) \to \mathbb{R}$ given by

$$
J(u) = \frac{\varepsilon}{2} \int_{\Omega} |\nabla u|^2 dx + \frac{1}{4\varepsilon} \int_{\Omega} (u^2-1)^2 dx.
$$

The extremals of $J$ satisfy the Euler-Lagrange equation (in the weak sense)

$$
-\varepsilon \Delta u + \frac{1}{\varepsilon} u(u^2-1) = 0 \qquad \qquad \text{in} \quad \Omega
$$

with the Neumann boundary condition

$$
\frac{\partial u}{\partial n} = 0 \qquad \qquad \text{on} \quad \partial \Omega.
$$

The **Allen-Cahn equation** is the $L^2$-gradient flow for $J$, which is governed by (here $u = u(t,x), t \geq 0, x \in \overline{\Omega}$) 

$$
u_t = \varepsilon \Delta u - \frac{1}{\varepsilon} u(u^2-1), \qquad \qquad \text{in} \quad \Omega
$$

$$
\frac{\partial u}{\partial n} = 0 \qquad \qquad \text{on} \quad \partial \Omega,
$$

$$
u(0,x) = u_0(x).
$$

By definition, the solution satisfies the initial condition $u(0,x) = u_0(x)$ and

$$
(u_t, v) = -\varepsilon(\nabla u, \nabla v) - \frac{1}{\varepsilon} (u(u^2-1),v)
$$

for all $v \in H^1(\Omega)$.

### Several important observations
- The following maximum principle holds for the Allen-Cahn equation.
If the initial condition satisifes $-1 \leq u_0 \leq 1$, then $-1 \leq u(t,x) \leq 1$ for all $t \geq 0$.
- The constant functions $-1$ and $+1$ are stable minimizers of the functional $J$.
- As $\varepsilon > 0$ is small, it is expected to observe *phase separation*. That is, after sufficient time, in large portions of the domain, the solution is close to either $+1$ or $-1$. The phases are separated by thin phase separation layers.

### Description of discretization procedure

Let $h>0$ and let $V_h \subset H^1(\Omega)$ be a finite dimensional subspace of $H^1(\Omega)$.
We use the following discretization algorithm: starting with $V_h \ni u^0_h \approx u_0$, find $(u^k_h \in V_h)_{k \geq 0}$ that satisfy

$$
(\frac{u^{k+1}_h-u^k_h}{\alpha_k}, v_h) = -\varepsilon(\nabla u^{k+1}_h, \nabla v_h) - \frac{1}{\varepsilon} (u^k_h((u^k_h)^2-1),v_h)
$$

for all $v_h \in V_h$. Here $\alpha_k > 0$ are the time discretization step sizes. The scheme can be rewritten as

$$
(u^{k+1}_h, v_h) + \alpha_k \varepsilon (\nabla u^{k+1}_h, \nabla v_h) = (u^k_h, v_h) - \frac{\alpha_k}{\varepsilon} (u^k_h((u^k_h)^2-1),v_h)
$$

for all $v_h \in V_h$. Some comments:

- We treat the diffusion term implicitly, and the nonlinear term explicitly.
- We use first-order $H^1$ finite element spaces.
- We use the backtracking line search and the Armijo rule (from Module 7) to choose the step sizes $\alpha_k$ in the iteration procedure. As the initial guess of the step size, we take a number of order $\varepsilon$.
- The stopping rule is

$$
\Vert \frac{u^{k+1}_h-u^k_h}{\alpha_k} \Vert = \Vert\nabla J(u^k_h)\Vert < tol,
$$

where $tol$ is a predetermined tolerance. The norm and $\nabla J$ above are taken in the sense of $L^2(\Omega)$.

### Results
We tested our algorithm on a variety of meshes. Input parameters: $\varepsilon = 0.4$ and exit tolerance `tol=0.05`.
The initial condition $u^0_h$ was taken to be a `GridFunction` with random coefficients generated from the uniform distribution on the interval (-0.95,0.95).

This repository contains pictures of initial data and the corresponding computed solutions at the exit times for various meshes. The repository also contains the meshes on which the tests were performed. The meshes are taken from MFEM.
