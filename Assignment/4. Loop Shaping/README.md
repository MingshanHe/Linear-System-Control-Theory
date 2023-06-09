<img src="LOGO\SNU.png" style="zoom:6%;" /><img src="LOGO\Seoul_national_university_logotype.svg.png" style="zoom:10%;" />



#### 1. The signal $u(t)$ is referred to as a *power* signal if its average power has a finite value. We define the function *pow* as the square root of the average power:

$$
pow(u)\coloneqq \begin{pmatrix}\lim_{T\to\infty}{\frac{1}{2T}\int^{T}_{-T}u(t)^2dt}\end{pmatrix}^{\frac{1}{2}}
$$

#### Although *pow* does not satisfy all the criteria of a norm, it can be used to establish an inclusion relation in the Venn diagram show in Fig.1a. Each set represents the set of signals that have a finite value according to the given norm definition (or *pow* definition).

<img src="Figure\1.PNG" style="zoom:80%;" />

##### (a) To prove the validity of Fig.1a, demonstrate the following statements:

###### 	i. If $\lVert{u}\rVert_2 <\infty$, then $u(t)$ is a *power signal* with $pow(u)=0$

​	Because of $\lVert{u}\rVert_2<\infty$, then we have $(\Sigma_{i=1}^n{u_i^2})^{\frac{1}{2}}<\infty$. And we can assume that term is a constant $C$ which is smaller than infinity. Then,
$$
pow(u)\coloneqq \begin{pmatrix}{\lim_{T\to\infty}\frac{1}{2T}\int^T_{-T}u(t)^2dt}\end{pmatrix}^{\frac{1}{2}}
= \begin{pmatrix}{\lim_{T\to\infty}(\frac{1}{2T})^{\frac{1}{2}}\int^T_{-T}\lVert{u}\rVert_2 dt}\end{pmatrix}
= \begin{pmatrix}{\lim_{T\to\infty}\frac{1}{2T}^{\frac{1}{2}}C}\end{pmatrix}
=0
$$

###### 	ii. If $\lVert{u}\rVert_1<\infty$ and $\lVert{u}\rVert_{\infty}<\infty$, then $\lVert{u}\rVert_2<\infty$

​	Because of $\lVert{u}\rVert_1<\infty$, we have $\min_{i\in R} u_i < \infty$. And $\lVert{u}\rVert_{\infty} < \infty$, we have $\max_{i\in R} u_i < \infty$. Then, we have,
$$
\lVert{u}\rVert_2 = (\Sigma_{i=1}^n{u_i^2})^{\frac{1}{2}}
$$
​	where the equation is bounded by the $L_1$ norm and $L_\infty$ norm,
$$
\infty\geq(\Sigma_{i=1}^n{u_i^2})^{\frac{1}{2}} \geq (\Sigma_{i=1}^n{\min u^2})^{\frac{1}{2}}\\
(\Sigma_{i=1}^n{u_i^2})^{\frac{1}{2}} \leq (\Sigma_{i=1}^n{\max u^2})^{\frac{1}{2}}<\infty
$$
​	Then we can easily know that $L_2$ norm is bounded, it has a finite value.

##### (b) Determine the location $(1,...,9)$ where $u$ should be included in Fig.1b. Given that $u(t) = 0$ for $t<0$, consider the following cases:

###### 	i. $u(t) = 1$ 

​		$pow(u) = 1$ for $\forall t$			5

###### 	ii. $u(t) = \left\{ \begin{aligned} \frac{1}{\sqrt{t}} && t\leq1 \\  0 && t >1 \end{aligned}\right.$

​		$pow(u) \leq 1 $						6

###### 	iii. $u(t) = \left\{ \begin{aligned} t^{-\frac{1}{4}} && t\leq1 \\  0 && t >1 \end{aligned}\right.$

​		$pow(u)\leq 1$						2

###### 	iv. $u(t)=\sum_{k=1}^\infty v_k(t)$ where $v_k(t)= \left\{ \begin{aligned} k && k<t<k+k^{-3} \\  0 && otherwise \end{aligned}\right.$

​		8

#### 2. Consider a plant $P(s) = 1/(s-1)$ with a unity feedback system. Suppose a disturbance $\omega(t) = asin(2t+\theta)$, with unknown amplitude $a$ and phase $\theta$, enters the plant as shown in Fig.2. Design a compensator of degree 3 that is proper (but not strictly proper) in such a way that the output asymptotically tracks any step reference input and rejects the disturbance. Place the poles at $-1\pm2j$ and $-2\pm1j$. Verify your results using MATLAB.

<img src="Figure\2.PNG" style="zoom:80%;" />

​	Assume that the disturbance $\omega(t)=0$, then we can compute the transfer function $G_1(s)$ with input signal $r(t)$ and output signal $y(t)$:
$$
G_1(s) = \frac{C(s)P(s)}{1+C(s)P(s)}
$$
​	And assume that the input signal $r(t) = 0$, then we can compute the transfer function $G_2(s)$ with the disturbance signal $\omega(t)$ and output signal $y(t)$:
$$
G_2(s) = \frac{P(s)}{1+C(s)P(s)}
$$
​	Then the output can be represented with above transfer functions:
$$
y(s) = G_1(s)r(s) + G_2(s)\omega(s)
$$
​	For the internal model, the $\phi(s) = \frac{1}{r(s)\omega(s)} = s^3+4s$, Then $B(s)/A(s)$ can be solved from:
$$
B(s)+A(s)(s-4)\phi(s) = F(s)
$$
​	where $F(s) = (s+(-1+2j))(s+(-1-2j))(s+(-2+1j))(s+(-2-1j))=s^4+6s^3+18s^2+30s+25$.

​	Because the compensator of degree is 3 and the degree of $\phi(s)$ is 3 now, thus $A(s) = 1$. And $B(s) = 10s^3+14s^2+46s+25$ by solving the equation. Then the compensator is:
$$
C(s) = \frac{B(s)}{A(s)\phi(s)}=\frac{10s^3+14s^2+46s+25}{s^3+4s}
$$
<img src="Figure\3.PNG" style="zoom:50%;" />