# SRQ2 Min - Minimization of Sum of Two Rayleigh Quotients

Dated 		06-20-2024

This folder contains the MATLAB codes used in the research paper:

>* Eigenvalue Backward Errors of Rosenbrock Systems and Optimization of Sums of Rayleigh Quotient*
by Shreemayee Bora, Ding Lu, Anshul Prajapati, and Punit Sharma.
Technical Report, 2024.


## Description

The eigenvalue backward error of a Rosenbrock matrix polynomial can be expressed as an SRQ2 minimization
```math
	\min_{x\in\mathbb C^n,\ \|x\|=1}
	f(x):= \frac{x^*A_1x}{x^*(\alpha_1 I+\beta_1 A_3)x} + \frac{x^*A_2x}{x^*(\alpha_2 I+\beta_2 A_3)x},
```
where $A_1,A_2,A_3\in\mathbb C^{n,n}$ are Hermitian and positive semidefinte matrices, $\alpha_i,\beta_i$ are given scalars satisfying $\alpha_i I+\beta_i A_3\succeq 0$, for $i=1,2$.

We adopt a nonlinear eigenvector approach for SRQ2 minimization, based on the fact that the solution $x_*$ of SRQ2 minimization can be characterized by a nonlinear eigenvector problem (NEPv) in the form of 
```math
 H(x) x = \lambda x,
```
where $H(x)$ is a Hermitian matrix-valued function and $\lambda$ is the smallest eigenvalue of $H(x)$. In our implementation, the NEPv is solved by *self-consistent-field iteration* (SCF) with adaptive level-shifting. This package contains all the example data and routines used in our study to facilitate its reproducibility.


## Contents

- Examples 
	- example1_localopt.m:	demonstrate local minimizers of SRQ2 and JNR minimization
	- example2_full.m:		demonstrate efficiency of NEPv approach against manifold optimization; full perturbation
	- example3_partial.m:	demonstrate efficiency of NEPv approach against manifold optimization; partial perturbation
- Private
	- buildSRQ2.m:	generate SRQ2 minimization for eigenvalue backward error of Rosenbrock systems
    - Solvers: 		runscf2.m (SCF iteration), runrtr.m (Riemannian trust-regions)
	- Other files: 	getinitial.m (generate initial vector), samplejnr3.m (sample boundary points of Joint Numerical Range), showjnr3.m (draw JNR), showcontour.m (draw contour of $g(y)=c$)


## External

The MATLAB toolbox Manopt is used for comparison. Please first download
and install the [Manopt](https://www.manopt.org/).


## Contact 

For questions, please contact Ding.Lu@uky.edu  
