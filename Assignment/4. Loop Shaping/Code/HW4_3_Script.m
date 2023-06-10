clc;
clear all;
close all;
syms a1 a2 b1 b2 b3
eqn1 = 2+a1+b1 == 7;
eqn2 = 3+2*a1+a2+b2 == 20;
eqn3 = 3*a1+2*a2+4-b1+b3==30;
eqn4 = 4*a1+3*a2-b2 == 24;
eqn5 = 4*a2-b3 == 8;
sol = solve(eqn1,eqn2, eqn3,eqn4, eqn5,a1,a2, b1, b2,b3);
sol.a1
sol.a2
sol.b1
sol.b2
sol.b3