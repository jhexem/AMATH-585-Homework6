y0 = 1;
f = @(t, y) y.^2 - sin(t) - cos(t).^2;

T = 8;
Nvals = [25, 50, 100, 200, 400, 800, 1600];

idx = 1;
N = Nvals(idx);

h = T / N;

[xFE, yFE] = forwardEuler(f, y0, h, N);
[xRK, yRK] = fourthOrderRK(f, y0, h, N);
%{
plot(xRK, yRK, "o")
hold on;
plot(xRK, cos(xFE), "-r", "LineWidth", 2)
title("Classical Fourth-Order Runge-Kutta Approximation vs True Solution")
legend("Classical Fourth-Order Runge-Kutta Approximation", "True Solution")
xlabel("t")
ylabel("y(t)")
%}

errorFE = zeros(1, 7);
errorRK = zeros(1, 7);
hvals = (T * ones(1, 7)) ./ Nvals;
ytrue = cos(T);

for j = 1:7
    N = Nvals(j);
    h = hvals(j);
    [xFE, yFE] = forwardEuler(f, y0, h, N);
    [xRK, yRK] = fourthOrderRK(f, y0, h, N);
    errorFE(j) = abs(yFE(end) - ytrue);
    errorRK(j) = abs(yRK(end) - ytrue);
end
%
loglog(hvals, errorFE, "-b", "LineWidth", 2)
hold on;
loglog(hvals, errorRK, "-r", "LineWidth", 2)
title("Errors for Forward Euler and Classical Fourth-Order Runge-Kutta")
legend("Forward Euler", "Runge-Kutta")
xlabel("h")
ylabel("error")
%}
function [x, y] = fourthOrderRK(f, y0, h, N)
    x = (h .* (0:N)).';
    y = [y0; zeros(N, 1)];
    for j = 2:(N+1)
        k1 = f(x(j-1), y(j-1));
        k2 = f(x(j-1) + (h/2), y(j-1) + (h/2).*k1);
        k3 = f(x(j-1) + (h/2), y(j-1) + (h/2).*k2);
        k4 = f(x(j-1) + h, y(j-1) + h.*k3);
        y(j) = y(j-1) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function [x, y] = forwardEuler(f, y0, h, N)
    x = (h .* (0:N)).';
    y = [y0; zeros(N, 1)];
    for j = 2:(N+1)
        y(j) = y(j-1) + h .* f(x(j-1), y(j-1));
    end
end