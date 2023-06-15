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
hold on
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