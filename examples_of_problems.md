The following are five typical examples of solving optimization numerical problems. For the solution code, see the path './util/'

1. Using the steepest gradient descent method to solve
   $$
   \min f(X)=x_1^2+25 x_2^2
   $$
   Initial value $X_0=[1,2]^T, \varepsilon=0.01$.

2. Using the Newton's method to solve
   $$
   \min f(X)=x_1^2+x_2^2-x_1 x_2-10 x_1-4 x_2+60
   $$
   Initial value $X_0=[0,1]^T, \varepsilon=0.01$.

3. Using the conjugate gradient method to solve
   $$
   \min f(X)=x_1^2+4 x_2^2
   $$
   Initial value $X_0=[1,1]^T, \varepsilon=0.01$.
4. Given the function 
   $$
   f(x)=x_1^2+x_2^2-1,
   $$
   at the current point $ x=[2,2]$, with a descent direction $ d=[-1,-1]^ T $, use the Armijo search method to obtain an acceptable step length ( $ \rho=0.1 $​ ).

5. Solve the unconstrained optimization problem using the trust region metho
   $$
   \min _{x \in R^2} f(x)=100\left(x_1^2-x_2\right)^2+\left(x_1-1\right)^2
   $$
   This problem has an exact solution $x^*=[1,1]^T, f\left(x^*\right)=0$.