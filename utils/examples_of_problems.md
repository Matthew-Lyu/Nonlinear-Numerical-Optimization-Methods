The following are five typical examples of solving optimization numerical problems. For the solution code, see the path './utils/'

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

6. Solving Using Least Squares Method. When studying the rate of a single-molecule chemical reaction, the following data were obtained:

   |      Index       | $\mathbf{1}$ | $\mathbf{2}$ | $\mathbf{3}$ | $\mathbf{4}$ | $\mathbf{5}$ | $\mathbf{6}$ | $\mathbf{7}$ | $\mathbf{8}$ |
   | :--------------: | :----------: | :----------: | :----------: | :----------: | :----------: | :----------: | :----------: | :----------: |
   | $\boldsymbol{t}$ |      3       |      6       |      9       |      12      |      15      |      18      |      21      |      24      |
   | $\boldsymbol{y}$ |     57.6     |     41.9     |     31.0     |     22.7     |     16.6     |     12.2     |     8.9      |     6.5      |

   Here, $t$ represents the time elapsed since the start of the experiment, and $y$ represents the amount of reactant produced at time $t$. Based on the above data, determine the empirical formula $y=f(t)$. 

7. Solve the Quadratic Programming Problem and Write the Lagrangian Function and Analysis Process
   $$
   \left\{
   \begin{array}{l}
   \min \quad x_1^2 + x_2^2 \\
   \text{s.t. } x_1 + x_2 = 1 \\
   x_2 \leq \frac{1}{4}
   \end{array}
   \right.
   $$

8. Solve the Following Problem Using the Exterior Penalty Function Method with Precision $10^{-8}$ and Initial Point $(10,10)$

   $$
   \begin{gathered}
   \min f(x) = x_1^2 + x_2^2 \\
   \text{s.t. } x_1 + x_2 = 2
   \end{gathered}
   $$

9. Solve the Following Problem Using the Interior Penalty Function Method with Precision $10^{-8}$ and Initial Point $(10,10)$

   $$
   \begin{gathered}
   \min f(x) = x_1^2 + x_2^2 \\
   \text{s.t. } x_1 - 1 \geq 0
   \end{gathered}
   $$

10. Use Genetic Algorithm to Find the Maximum Value of the Following Function in the Specified Constraint Set:

    $$
     f(x) = x \sin (10x) + 1
    $$

     The constraint set is $\{x \in \mathbb{R} : x \in [-1, 2]\}$. The solution should be accurate to 6 decimal places.