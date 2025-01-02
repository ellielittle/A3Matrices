# A3Matrices

NOTE! This code requires functions from the folding table code (link TBD). After downloading this code download the folding table code and copy the file "ft_main.m" into the folder with this matrix code so that it is loaded when the code is run.

This code calculates the matrices of the matrice representations $\pi_{J,\mathsf{v}}$ described in Section 5 of https://arxiv.org/abs/2212.10781. These matrices are used in Example 7.13 of https://arxiv.org/abs/2406.07004 for the case when $J = \\{1,2\\}$.

The input is an element of the extending affine group of type $A_3$, notated $\widetilde{W}$. Let $\sigma$ be such that $\sigma(i) = i+1 \mod 4$, then $\widetilde{W} = W\rtimes \Sigma$ where $W$ is the affine coxeter group of type $A_3$ and $\Sigma = \{e,\sigma,\sigma^2,\sigma^3\}. Thus, elements of $\widetilde{W}$ are of the form $s_{i_1}\dots s_{i_l}\sigma^k$ for some $0\leq k\leq 3$ and $i_1,\dots, i_l\in I\cup\\{0\\}$ (where $I\cup\\{0\\} is the index set of the generators of $W$). 

To print the matrix, and the leading matrix (see Definition 6.3 of https://arxiv.org/abs/2212.10781), of an element of $\widetilde{W}$ the command is printmatrixandleading(t,J) where $J$ is a set and $t$ is of the form $<[i_{1},\dots, i_l], k>$ (corresponding to an element of the form $s_{i_1}\dots s_{i_l}\sigma^k$). 

The $z$ indeterminants correspond to the $z$ indeterminants defined in Section 2.5 of https://arxiv.org/abs/2406.07004, so are determined by the partition corresponding to $J$. Note that the matrices the code creates does not factor in the relations between the $z$ indeterminants so will not cancel terms nicely (that is no ideal is created to quotient). 
