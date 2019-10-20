%%ODE solver for enzyme dynamics (iGEM, bioreactor)
clear all; clc;
%Initial concentration of PET
so = 100; %mM of PET
%%Initial conditions
y0 = [so, 0, 0, 0, 0]; %[s(0), c(0), p_1(0), c'(0), p_2(0)]
%%Time interval for solution
%t = linspace(0, 30, 100);

%%Solution
[t, y] = ode45(@f, [0, 30] , y0);
plot(t, y(:, 1))
hold on
plot(t, y(:, 2))
hold on
plot(t, y(:, 3))
hold on
plot(t, y(:, 4))
hold on
plot(t, y(:, 5))
legend('PET', 'PET + PETase' , 'MHET', 'MHET + MHETase complex', 'EG + TPA')
xlabel('Time (A.U)')
ylabel('Concentration (A.U)')
title('Concentration vs Time (Tank 1)')
%%ODEs
function model = f(t,y)

S = y(1);
C1 = y(2);
P1 = y(3);
C2 = y(4);
P2 = y(5);

k_1 = 2;   
E0_1 = 5; %mM 
k_n1 = 1;  
k_2 = 2;   
k_3 = 2;   
k_n3 = 1;  
E0_2 = E0_1;  %mM 
k_4 = 1;    

dSdt = -k_1*(E0_1-C1)*S + k_n1*C1;
dC1dt = k_1*(E0_1-C1)*S - (k_n1+k_2)*C1;
dP1dt = k_2*C1 - k_3*P1*(E0_2-C2) + k_n3*C2;
dC2dt = k_3*(E0_2-C2)*P1 - (k_n3+k_4)*C2;
dP2dt = k_4*C2;

model = [dSdt; dC1dt; dP1dt; dC2dt; dP2dt];
end


