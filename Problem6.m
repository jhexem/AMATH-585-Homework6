y0 = [20; 20];
tspan = [0; 20];

[T, Y] = ode45(@lotkaVolterra, tspan, y0);
%{
plot(T, Y(:, 1), '-r', "LineWidth", 2)
hold on;
plot(T, Y(:, 2), '-b', "LineWidth", 2)
title("Lotka-Volterra Predator-Prey Equations")
legend("R(t)", "F(t)")
xlabel("t")
ylabel("Population")
%}
%
plot(Y(:, 1), Y(:, 2), '-r', "LineWidth", 2)
title("Lotka-Volterra Predator-Prey Equations")
xlabel("R(t)")
ylabel("F(t)")
%}
function dydt = lotkaVolterra(t, y)
    dydt = zeros(2, 1);
    dydt(1) = (1 - 0.02*y(2))*y(1);
    dydt(2) = (-1 + 0.03*y(1))*y(2);
end