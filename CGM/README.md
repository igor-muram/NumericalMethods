# CGM
Solving linear systems of equations using the conjugate gradient method and using a locally optimal scheme.<br>
Solution of systems of linear equations without preconditioning of the coefficient matrix.<br>
Implemented LU factorization.<br>
Hilbert matrix generator of any dimension in sparse format.<br>
Implemented three types of diagonal preconditioning for the coefficient matrix for CGM and LOS.<br>
Implemented LU preconditioning of the coefficient matrix for CGM and LOS.<br>
The matrix of coefficients is stored in a sparse row-column format.<br>
The result is the vector x (vector of unknowns), calculated by one of two methods, the number of iterations for the selected iterative method, and the iterative residual.<br>
The selection of the preconditioning of the coefficient matrix is ​​carried out in the kuslau file.<br>