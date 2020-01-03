# NewtonsMethod
The solution of the system of equations by the Newton's method (tangent method).<br>
The matrix of coefficients is stored in a sparse row-column format.<br>
If the number of equations and the number of variables do not match, then the equation that has the smallest absolute value is deleted.<br>
Automatic selection of the relaxation parameter for guaranteed convergence to the solution.<br>
It is possible to use the calculated results to visually represent the convergence process of the Newton's method. To do this, after computing the solution, you need to run the visualizer, which is located in this repository.<br>
The result is a point that is a solution to a system of equations.