<img src="LOGO\SNU.png" style="zoom:6%;" /><img src="LOGO\Seoul_national_university_logotype.svg.png" style="zoom:10%;" />

### 1. Consider a CT-LTI System:

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

​	For the open-loop observer, this means that the observer using the same state space equation system to tracking the original system. Because there has some situations that the states of the system are very expansive or hard to measure.                                                                                                  

