%Tank 3, bioreactor

%%ODE solver for enzyme dynamics (iGEM, bioreactor)

%Initial concentration of PET
so = 5;

%%Initial conditions
y0 = [0, 0, 0]; %[s(0), c(0), p_1(0), c'(0), p_2(0)]

%%Time interval for solution
time = linspace(1,10, 100);

%%Solution
[t, y] = ode45(@f, time, y0);

%plot
plot(time, y(:, 1), '-s', time, y(:, 2), '-o', time, y(:, 3), '*')
legend('NH4OH', 'TPA', 'Diammonium Salt')
xlabel('Time')
ylabel('Concentration')
title('V = 15, NH4OH flux = 5, Tank 2 flux = 5, reaction rate = 15')


%%ODEs
function model = f(t,y)
%Variables in Andre's notes: y(1) = s, y(2) = c, y(3) = p_1, y(4) = c',
%y(5) = p_

V = 15;
B_in = 5;
L_in = 5;
k = 15;

model = zeros(3,1);
model(1) = B_in/V - k*y(2)*y(1);
model(2) = L_in/V - k*y(2)*y(1);
model(3) = k*y(2)*y(1);
end