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

#### 3. Given the plant transfer function $P(s)$, implement the model $H_0$:

$$
P(s) = \frac{s^2-1}{s^3+2s^2+3s+4}, \space H_0(s) = \frac{(s-1)(2s+1)}{(s+2)^2(s^2+2s+2)}
$$

#### by designing a feedforward pre-compensator $C_1(s)$ and a feedback controller $C_2(s)$. Determine if the resulting system is stable and check desired model matching is achieved by using MATLAB.

<img src="Figure\4.PNG" style="zoom:70%;" />

​	$H_0(s)$ is implementable:
$$
\frac{H_0}{N_P}=\frac{2s+1}{(s+1)(s+2)^2(s^2+2s+2)}=\frac{\bar{E}}{\bar{F}}
$$
​	Then we need to find $M_c$ and $N_2$ that satisfy:
$$
M_pM_c+N_pN_2=\bar{F}
$$
​	Using `solve` function in MATLAB, we can solve this equation with,

```matlab
syms a1 a2 b1 b2 b3
eqn1 = 2+a1+b1 == 7;
eqn2 = 3+2*a1+a2+b2 == 20;
eqn3 = 3*a1+2*a2+4-b1+b3==30;
eqn4 = 4*a1+3*a2-b2 == 24;
eqn5 = 4*a2-b3 == 8;
sol = solve(eqn1,eqn2, eqn3,eqn4, eqn5,a1,a2, b1, b2,b3);
```

​	Then, $a_1 = 4.5$, $a_2 = 3.5$, $b_1 = 0.5$, $b_2 = 4.5$, $b_3 = 6$.

​	

<img src="Figure\5.PNG" style="zoom:50%;" />



#### 4. Consider the pitch rate control of aircraft $P(s)$ where reference $r$ is pitch rate command, and output $y$ is pitch rate of the aircraft, $P_b(s)$ is bending mode is considered model:

$$
P(s) = \frac{s+1}{s^2+7s+25}\\
P_b(s) = P(s)\frac{s^2+3s+30^2}{s^2+0.9s+45^2}
$$

##### (a)	Consider a controller $C(s)$ obtained using the loop transfer function $L(s) = \omega_c/s$ without considering bending. Use the maximum value of $\omega_c$ that satisfies $\vert L_b(j\omega)\vert < 0.5$ for all $\omega \geq 45$. Plot the loop shape $\vert L_b(j\omega)\vert$ considering bending using MATLAB when using this controller.

​	I plotted three systems:

		1. without lead-lag compensator and just original $P(s)$
		1. with lead-lag compensator and just original $P(s)$
		1. with lead-lag compensator and considered model $P_b(s)$

​	And I also plotted the Bode figure in the MATLAB. I tried to design the parameter $w_c$ which is desired the condition in this question.

```matlab
clear
clc
close all
s = tf('s');
P = (s+1)/(s^2+7*s+25);
Pb = P * (s^2+3*s+900)/(s^2+0.9*s+2025);
%% Lead-lag controller
wc = 0.5;
F = wc/s;
%% plot
figure(1)
margin(P)
hold on
margin(P*F)
hold on
margin(Pb*F)
legend('w/o lead-lag controller','w lead-lag controller','w lead-lag controller w bending mode');
hold off

figure(2)
stepplot(feedback(P,1))
stepinfo(feedback(P,1))
hold on
stepplot(feedback(P*F,1))
stepinfo(feedback(P*F,1))
hold on
stepplot(feedback(Pb*F,1))
stepinfo(feedback(Pb*F,1))
legend('w/o lead-lag controller','w lead-lag controller','w lead-lag controller w bending mode');
```

<img src="Figure\6.png" style="zoom:150%;" />

<img src="Figure\7.png" style="zoom:150%;" />

##### (b) Design a proper controller $C_b(s)$ to achieve a loop shape with the following properties: a) The loop shape should have a similar value to the original loop shape $\omega/c$ at low frequencies, and b) The loop shape should satisfy $\vert L_b(j\omega)\vert<0.5$ for all $\omega \geq 45$. Find the cut off frequency $\omega_c$ for controller $C_b$ and compare the loop shape obtained in  Problem (4a) using MATLAB.

​	I am so sorry for this question. I have no idea.







​																																																									<img src="C:/Users/BRL/Documents/GitHub/Control-System-1/Assignment/3. State Feedback And Kalman Filter/LOGO/signature.png" style="zoom:120%;" />

​																																																										淡泊名利，矢志不渝；

​																																																										上下求索，苦心孤诣；

​																																																										百折不挠，夜以继日；

​																																																										宁静致远，知行合一；