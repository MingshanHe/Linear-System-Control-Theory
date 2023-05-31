<img src="LOGO\SNU.png" style="zoom:6%;" /><img src="LOGO\Seoul_national_university_logotype.svg.png" style="zoom:10%;" />

### **1. Consider a CT-LTI System:**

$$
\dot{x} = \begin{bmatrix}1&-2& 1\\1&1&3\\-1&4&0\end{bmatrix}x+\begin{bmatrix}1\\-1\\0\end{bmatrix}u
$$

$$
y = \begin{bmatrix}1&0&2\end{bmatrix}x
$$

#### (1) Design state feedback controller and state observer

##### 	a) Design the state feedback system that have $-5, -5+3j, -5-3j$ as its eigenvalues.

​	To design a state feedback system with desired eigenvalues, two distinct approaches are commonly employed: one involving the utilization of the "place" function and another without it. The process begins by establishing the state space equation, enabling the computation of eigenvalues and the characteristic polynomial equation for the original system.

​	By formulating the system's dynamics using matrices A, B, C, and D, the eigenvalues can be determined from the A matrix, providing valuable insights into the system's stability characteristics. Additionally, the characteristic polynomial equation is derived from these eigenvalues, offering a mathematical representation of the system's behavior.

```matlab
%Define State Space Equation System
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0; 
sys = ss(A, B, C, D);
eigenvals = eig(A);
char_poly = poly(eigenvals);
```

​	In this section, the computed eigenvalues of the system are identified as $-2.0791$, $0.1215$, and $3.9576$, leading to a characteristic polynomial equation of $1*s^3 - 2*s^2 - 8*s + 1$. To design the appropriate gain $K$ for the state feedback controller, it is necessary to compute the transformation matrix $P$, as mentioned in the previous section. Additionally, the transformed gain $\bar{K}$ is computed to transition from the original eigenvalues to the desired eigenvalues. Finally, the final gain $K$ can be calculated using $\bar{K}$ and $P$. This section code:

```Matlab
% Compute State Feedback Gain K w/o 'place Funtion
C_ = [B A*B A*A*B];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
K_bar = [17 92 169];
K = K_bar*P;
```

​	For compute state feedback gain $K$, MATLAB provide the function `place` to compute:

```matlab
% Compute State Feedback Gain K w 'place Funtion
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);
```

​	After performing the computational analysis, the determined state feedback gain, denoted as $K$, is found to be $-17.5570$, $-34.5570$, and $-35.7342$. This selected gain value effectively transforms the eigenvalues of the original system to the desired values, thereby achieving the desired system stability.  The following result presents a comparative analysis of the step response between the original system and the state feedback controller system. The eigenvalues of the original system exhibit a positive portion, indicating instability. Conversely, the state feedback controller system successfully achieves stability by ensuring all eigenvalues reside in the negative domain, as evidenced by its stable step response.

<img src="Figure\1.png" style="zoom:100%;" />

##### b) Design three state estimators: open-loop observer, two Luenberger observers that each have following sets of eigen values, $-1, -2, -3$ and $-0.01, -0.02, -0.03$. Plot the output responses from the initial condition $x(0) = \begin{bmatrix}1&1&1\end{bmatrix}^T, \hat{x}(0)=\begin{bmatrix}0&0&0\end{bmatrix}^T$ and discuss about them.

​	For the **open-loop observer**, this means that the observer using the same state space equation system to tracking the original system. Because there has some situations that the states of the system are very expansive or hard to measure. In textbook, it mentioned that there are two disadvantages in using an open-loop estimator. First, the initial state must be computed and set each time we use the estimator. This is very inconvenient and may be expansive to measure. Second, and more seriously, if the matrix A has eigenvalues with positive real parts (the system is unstable) then even for a very small difference between $x(t_0)$ and $\hat{x}(t_0)$ for some $t_0$, which may be caused by disturbance or imperfect estimation of the initial state, the difference between $x(t)$ and $\hat{x}(t)$ will grow with time.                    

```matlab
% Mainly for the State Feedback and Track Controller
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_estimate0 = [0; 0; 0];% Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(3, length(t));

x_estimate = zeros(3, length(t));
y_estimate = zeros(3, length(t));
x(:, 1) = x0;
x_estimate(:, 1) = x_estimate0;
for i = 2:length(t)
    % Update the state using the system equations
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    % Update the estimated state using the observer equations
    x_hat(:, i) = x_hat(:, i-1) + dt * (A*x_hat(:, i-1) + B*u);
end                                                           
```

​	As the result of this section code, the figure has shown three states of the original system and observer with different initial state. Because the matrix $A$ of the original system has the real parts, and the small different from the initial states. The observer can not detect real states and track it successfully.

<img src="Figure\2.png" style="zoom:100%;" />

​	For the **closed-loop observer**, we shall modify the estimator making the output $y(t)=\pmb{cx}(t)$ compared with $\pmb{c\hat{x}}(t)$. And the difference passing through an $n\times1$ constant gain vector $L$. And I compute the gain $L$ with the similar method of gain $K$. As the condition, the eigen values set by the closed-loop observer.

```matlab
% Compute Observer Gain K w/o 'place Funtion
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Compute State Feedback Gain L w/o place Function
C_ = [C' A'*C' A'*A'*C'];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
L_bar = [8 19 5];
L = (L_bar*P)';
% Compute State Feedback Gain L and K w 'place Funtion
desired_poles = [-0.01 -0.02 -0.03];
L_placed = place(A', C', desired_poles)';
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);
```

​	After check the $L$ and $L_{placed}$, I use the gain and design the closed-loop observer. Then check the output response of different conditions. 

```matlab
% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_estimate0 = [0; 0; 0];     % Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));
x_estimate = zeros(3, length(t));
y_estimate = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
x_estimate(:, 1) = x_hat0;
y_estimate(:, 1) = C*x_estimate0;

for i = 2:length(t)
    u = -K_placed * x_hat(:, i-1);
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    y(i) = C*x(:, i);
    x_estimate(:, i) = x_hat(:, i-1) + dt * ((A - L*C)*x_hat(:, i-1) + B*u + L*y(i));
    y_estimate(i) = C*x_estimate(:, i);
end
```

​	Then the result has shown that the closed-loop observer can track the real states successfully and the output response is correct.

<img src="Figure\3.png" style="zoom:100%;" />

#### (2) Design the reference tracking controller so that $\lim_{t\to\infty}y(t)=1$. Plot the output response from the initial condition $x(0) = \begin{bmatrix}1&1&1\end{bmatrix}^T$.

​	Consider the reference tracking problem: $y(t)\to r(t)$.  For compute the satisfy conditions:
$$
0 = (AN_x+BN_u)r\\
r = (CN_x+DN_u)r
$$
​	And for the computation,
$$
N_x = \begin{bmatrix}8&2&-3.5\end{bmatrix}; N_u = -0.5
$$
​	Then for the reference controller, we can incorporate feedforward action into the feedback control, using the input
$$
u = -K*x+(K*N_x+N_u)*r
$$
​	From the result, the figure has shown the output response can successfully track the input reference signal.

```matlab
dt = 0.01;              % Time step
t_end = 10;   
K_Ref = [17 35 35];
R = 1;

t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
for i = 2:length(t)
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*(K_Ref*x(:, i-1)-85*R));
    y(:, i) = C*x(:, i); 
end
```

<img src="Figure\4.png" style="zoom:100%;" />

#### (3) Design the state feedback controller and state observer for the system, and plot $x(t)$, $\hat{x}(t)$.

​	Consider sinusoidal disturbance $d = \xi sin(\omega t)$ with unknown magnitude $\xi$ and frequency $\omega$:
$$
\dot{x}=Ax + B(u+d)\\
y = Cx
$$
​	where $d$ is a disturbance, only know to satisfy differential equation noise model:
$$
\dot{\xi}=A_d\xi, \space d=C_d\xi
$$
​	And for the sinusoidal signal,
$$
A_d = \begin{bmatrix}0&-\omega\\ \omega & 0\end{bmatrix},C_d=\begin{bmatrix}1&0\end{bmatrix}
$$
​	Thus the total(augmented) system dynamics is,
$$
{\begin{pmatrix}\dot{x}\\ \dot{\xi}\end{pmatrix} = \begin{bmatrix}A&BC_d\\ 0&A_d\end{bmatrix}\begin{pmatrix}{x}\\ {\xi}\end{pmatrix}+\begin{bmatrix}B\\ 0\end{bmatrix}u}
\\ y = \begin{bmatrix}C&0\end{bmatrix}\begin{pmatrix}{x}\\ {\xi}\end{pmatrix}
$$
​	I design a state observer for this system is,
$$
{\begin{pmatrix}\dot{\hat{x}}\\ \dot{\hat{\xi}}\end{pmatrix} = \begin{bmatrix}A&BC_d\\ 0&A_d\end{bmatrix}\begin{pmatrix}\hat{x}\\ \hat{\xi}\end{pmatrix}+\begin{bmatrix}B\\ 0\end{bmatrix}u+\begin{bmatrix}L_1\\L_2\end{bmatrix}(y-C\hat{x})}
$$
​	The code is,

```matlab
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
Ad = [0, -2; 2, 0];
Cd = [1, 0];
A_aug = [A B*Cd; zeros(2,3) Ad];
B_aug = [B; zeros(2,1)];
C_aug = [C zeros(1,2)];
% Compute State Feedback Gain L and K w 'place Funtion
desired_poles = [-4 -4+0.5i -4-0.5i -4+2i -4-2i];
L_placed = place(A_aug', C_aug', desired_poles)';
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

R = 1;
dt = 0.01
t = 0:dt:10;
x = zeros(5, length(t));
x_estimate = zeros(5, length(t));
y = zeros(1, length(t));
y_estimate = zeros(1, length(t));

% Define the initial conditions
x(:, 1) = [1; 1; 1; 5; 5];
x_estimate(:,1) = [0; 0; 0; 5; 5];
y(:, 1) = C_aug*x(:, 1) ; 
y_estimate(:,1) = C_aug*x_estimate(:,1);

for i = 2:length(t)
    u = -K_placed * x_estimate(1:3, i-1);
    x(:, i) = x(:, i-1) + dt * (A_aug*x(:,i-1)+ B_aug*u);
    y(:, i) = C_aug*x(:, i); 
    x_estimate(:, i) = x_estimate(:, i-1) + dt * (A_aug*x_estimate(:, i-1) + B_aug*u + L_placed*(y(i-1)-C_aug*x_estimate(:,i-1)));
end
```

<img src="Figure\5.png" style="zoom:100%;" />

<img src="Figure\6.png" style="zoom:100%;" />

### **2. Consider a following cart-pole system with mass of the cart and the particle on the massless pole $m_c$ and $m_p$ respectively.**

<img src="Figure\7.png" style="zoom:50%;" />

​	The dynamics of the system can be written as:
$$
H(q)\ddot{q}+C(q,\dot{q})\dot{q}+G(q)=Bu
$$
​	where $q:=\begin{bmatrix}x&\theta\end{bmatrix}^T$, and the above matrices are:
$$
H(q) =\begin{bmatrix}m_c+m_p&m_plcos(\theta)\\m_plcos(\theta)&m_pl^2\end{bmatrix},\\
C(q,\dot(q)) =\begin{bmatrix}0&-m_pl\dot{\theta}sin(\theta)\\ 0&0\end{bmatrix},\\
G(q) = \begin{bmatrix}0\\m_pglsin(\theta)\end{bmatrix},\\
B = \begin{bmatrix}1\\0\end{bmatrix}
$$
​	Let $m_c = m_p = l = 1$, $g = 10$ for gravity. Then the linearized dynamics can be achieved around the fixed point of interest with state $\begin{bmatrix}q&\dot{q}\end{bmatrix}^T$:
$$
\bar{A}=\left.{\begin{bmatrix}0&I\\-H^{-1}\frac{\partial G}{\partial q} & -H^{-1}C\end{bmatrix}}\right|_{x=\bar{x},u=\bar{u}},
\bar{B}=\left.{\begin{bmatrix}0\\H^{-1}B\end{bmatrix}}\right|_{x=\bar{x},u=\bar{u}}
$$

#### 	Design a LQR controller to stabilize the system around the fixed point, and show your result using MATLAB from the initial condition $q=\begin{bmatrix}0& 0.01\end{bmatrix}^T$, $\dot{q}=\begin{bmatrix}0& 0\end{bmatrix}^T$ with $Q=I$, $R=1$.
