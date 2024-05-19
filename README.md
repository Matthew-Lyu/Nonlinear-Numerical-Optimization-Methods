# Optimization-Theory-and-Methods

本仓库是对最优化方法的介绍和MATLAB代码复现。

### 1. 最速梯度下降法

**步0** 选取初始点 $x_0 \in \mathbb{R}^n$, 设置误差 $0 \leq \varepsilon \ll 1$. 令迭代次数 $k:=0$.

**步1** 计算梯度 $g_k=\nabla f\left(x_k\right)$. 若 $\left\|g_k\right\| \leq \varepsilon$, 算法停止, 输出 $x_k$ 作为近似最优解.

**步2** 取方向为负梯度方向 $d_k=-g_k$.

**步3** 由线搜索技术(我选择的算法是Armijo算法)确定步长因子 $\alpha_k$.

**步4** 令 $x_{k+1}:=x_k+\alpha_k d_k, k:=k+1$​, 转步 1 .

### 2. 阻尼牛顿法

**步0** 给定误差阈值 $\varepsilon$、阻尼因子 $\delta$ 。初始化点 $x_0 \in \mathbb{R}^n$，设置迭代次数 $k：=0$。

**步1** 计算当前点的梯度 $g_k=\nabla f(x_k)$。如果梯度的范数 $\|g_k\|$ 小于等于预设的误差阈值 $\varepsilon$，则算法终止，并输出近似最优解 $x^* \approx x_k$。

**步2** 计算当前点的海森矩阵 $G_k=\nabla^2 f(x_k)$，并求解线性方程组 $G_k d = -g_k$，得到方向 $d_k$。

**步3** 寻找满足Armijo条件的步长。定义 $m_k$ 为满足以下不等式的最小非负整数 $m$：

```math
f\left(x_k+\delta^m d_k\right) \leq f\left(x_k\right)+\sigma \delta^m g_k^T d_k
```

**步4**更新下一个点为 $x_{k+1} = x_k + \alpha_k d_k$，增加迭代次数 $k := k + 1$，然后返回步 1。

### 3. 修正牛顿法

**步0** 给定参数 $\delta \in(0,1), \tau \in[0,1], \sigma \in(0,0.5)$, 终止误差 $0 \leq \varepsilon \ll 1$. 初始点 $x_0 \in \mathbb{R}^n$. 令 $k:=0$.

**步1** 计算 $g_k=\nabla f\left(x_k\right), \mu_k=\left\|g_k\right\|^{1+\tau}$. 若 $\left\|g_k\right\| \leq \varepsilon$, 算法终止, 输出 $x_k$ 作为近似极小点.

**步2** 计算海森矩阵 $G_k=\nabla^2 f\left(x_k\right)$. 解线性方程组

```math
\left(G_k+\mu_k I\right) d=-g_k
```

得解 $d_k$.

**步3** 令 $m_k$ 是满足下列不等式的最小非负整数 $m$ :

```math
f\left(x_k+\delta^m d_k\right) \leq f\left(x_k\right)+\sigma \delta^m g_k^T d_k .
```

令 $\alpha_k=\delta^{m_k}, x_{k+1}:=x_k+\alpha_k d_k$.

**步4** 令 $k:=k+1$, 转步 1 .

### 4. 共轭梯度算法

**步0** 给定误差阈值 $0 \leq \varepsilon \ll 1$ 和初始点 $x_0$. 计算 $g_0=\nabla f\left(x_0\right)$. 令 $k:=0$.

**步1** 若 $\left\|g_k\right\| \leq \varepsilon$, 停算, 输出 $x^* \approx x_k$.

**步2** 计算搜索方向 $d_k:$

```math
d_k= \begin{cases}-g_k, & k=0, \\ -g_k+\beta_{k-1} d_{k-1}, & k \geq 1,\end{cases}
```

其中当 $k \geq 1$ 时, 其中 $\beta_{k-1}=\frac{g_k^T g_k}{g_{k-1}^T g_{k-1}}$（该式由 Fletcher 和 Reeves 给出的, 故称之为 Fletcher-Reeves 公 式, 简称 FR 公式）

**步3** 利用精确线搜索方法确定搜索步长 $\alpha_k$.

**步4** 令 $x_{k+1}:=x_k+\alpha_k d_k$, 并计算 $g_{k+1}=\nabla f\left(x_{k+1}\right)$.

**步5** 令 $k:=k+1$​​, 转步 1 .

### 5. Armijo准则

**步0** 给定 $\beta \in(0,1), \sigma \in(0,0.5)$. 令 $m:=0$.

**步1** 若不等式

```math
f\left(x_k+\beta^m d_k\right) \leq f\left(x_k\right)+\sigma \beta^m g_k^T d_k
```

成立, 置 $m_k:=m, x_{k+1}:=x_k+\beta^{m_k} d_k$, 算法终止. 否则, 转步 2 .

**步2** 令 $m:=m+1$, 转步 1 .

### 6. 信赖域方法

用信赖域方法求解无约束优化问题

```math
\min _{x \in \mathbb{R}^n} f(x)
```

设 $x_k$ 是第 $k$ 次迭代点. 记 $f_k=f\left(x_k\right), g_k=\nabla f\left(x_k\right), B_k$ 是 Hesse 阵 $\nabla^2 f\left(x_k\right)$ 的第 $k$ 次近似, 则第 $k$ 次迭代步的信赖域子问题具有如下形式:

```math
\begin{array}{ll}
\min & q_k(d)=g_k^T d+\frac{1}{2} d^T B_k d, \\
\text { s.t. } & \|d\| \leq \Delta_k
\end{array}
```

其中 $\Delta_k$ 是信赖域半径, $\|\cdot\|$ 是任一种向量范数, 通常取 2 -范数或 $\infty$-范数. 设子问题的最优解为 $d_k$, 定义 $\Delta f_k$ 为 $f$ 在第 $k$ 步的实际下降量:

```math
\Delta f_k=f_k-f\left(x_k+d_k\right)
```

$\Delta q_k$ 为对应的预测下降量:

```math
\Delta q_k=q_k(0)-q_k\left(d_k\right) .
```

再定义它们的比值为

```math
r_k=\frac{\Delta f_k}{\Delta q_k} .
```

一般地, 我们有 $\Delta q_k>0$. 因此, 若 $r_k<0$, 则 $\Delta f_k<0, x_k+d_k$ 不能作为下一个迭代点, 需要缩小信赖域半径重新求解子问题. 若 $r_k$ 比较接近 1 , 说明二次模型与目标函数在信赖域范围内有很好的近似, 此时 $x_{k+1}:=x_k+d_k$ 可以作为新的迭代点, 同时下一次迭代时可以增大信赖域半径. 对于其他情况, 信赖域半径可以保持不变. 下面给出求解无约束优化问题信赖域方法的一般框架.

**信赖域方法的一般框架如下：**

**步0** 选取初始参数 $0 \leq \eta_1<\eta_2<1,0<\tau_1<1<\tau_2, 0 \leq \varepsilon \ll 1 . x_0 \in \mathbb{R}^n$.取定 $\tilde{\Delta}>0$ 为信赖域半径的上限, 初始信赖域半径 $\Delta_0 \in(0, \tilde{\Delta}]$. 令 $k:=0$.

**步1** 计算 $g_k=\nabla f\left(x_k\right)$. 若 $\left\|g_k\right\| \leq \varepsilon$, 停止迭代.

**步2** 求解子问题 $(11)$ 的解 $d_k$.

**步3** 按 $(14)$ 式计算 $r_k$ 的值.

**步4** 校正信赖域半径.

```math
\Delta_{k+1}:=
\begin{cases}
\tau_1 \Delta_k, & \text { 若 } r_k \leq \eta_1, \\
\Delta_k, & \text { 若 } \eta_1<r_k<\eta_2, \\
\min \left\{\tau_2 \Delta_k, \tilde{\Delta}\right\}, & \text { 若 } r_k \geq \eta_2,\left\|d_k\right\|=\Delta_k .
\end{cases}
```

**步5** 若 $r_k>\eta_1$, 则令 $x_{k+1}:=x_k+d_k$, 更新矩阵 $B_k$ 到 $B_{k+1}$, 令 $k:=k+1$,转步 1 . 否则 $x_{k+1}:=x_k$, 令 $k:=k+1$, 转步 2 .

**对于子问题的求解我使用了光滑牛顿法，算法步骤如下：**

**步0** 选取 $\delta, \sigma \in(0,1), \mu_0>0, \lambda_0 \geq 0 . d_0 \in R^n$, 置 $z_0=\left(\mu_0, \lambda_0, d_0\right)$, $\bar{z}=\left(\mu_0, 0,0\right)$. 选取 $\gamma \in(0,1)$ 使 $\gamma \mu_0<1$ 及 $\gamma\left\|H\left(z_0\right)\right\|<1$. 令 $j:=0$.

**步1** 如果 $\left\|H\left(z_j\right)\right\|=0$, 算法终止; 否则, 计算 $\beta_j=\beta\left(z_j\right)$.

**步2** 求解下列方程组得解 $\Delta z_j=\left(\Delta \mu_j, \Delta \lambda_j, \Delta d_j\right)$,

```math
H\left(z_j\right)+H^{\prime}\left(z_j\right) \Delta z_j=\beta_j \bar{z} .
```

**步3** 设 $m_j$ 为满足下式的最小非负整数:

```math
\left\|H\left(z_j+\delta^{m_j} \Delta z_j\right)\right\| \leq\left[1-\sigma\left(1-\beta \mu_0\right) \delta^{m_j}\right]\left\|H\left(z_j\right)\right\| .
```

令 $\alpha_j:=\delta^{m_j}, z_{j+1}=z_j+\alpha_j \Delta z_j \text {. }$

**步4** 令 $j:=j+1$, 转步 1 .
