% clear all variables, leave not witnesses
clear('all')

% initialize system variables
Pe = 5.0;       % Peclet number; do not use values > 1000
Da = 5.0;       % Dalton number
N = 20;         % number of control volumes

% build matrices and vectors
A = zeros(N,N);
b = zeros(N,1);

% set some auxiliary variables
dz = 1.0 / N;
beta = 1.0 / (Pe * dz);
Sp = Da * dz;
Su = -1.0;

% populate matrix
for i = 2:N-1
    A(i,i-1) = beta + 0.5;
    A(i,i+1) = beta - 0.5;
    A(i,i) = -A(i,i-1) - A(i,i+1) - Sp;
end

% set boundary values
A(1,2) = beta - 0.5;
A(1,1) = -A(1,2) - Sp + Su;
A(N,N-1) = beta + 0.5;
A(N,N) = -A(N,N-1) - Sp;

% set vector solution
b(1) = Su;

% solve matrix-vector equation
u = A\b;

% plot result
plot(linspace(0.5 * dz, 1.0 - 0.5 * dz, N), u, 'o')
hold on;

% plot analytical solution
z = linspace(0, 1, 100);
beta = sqrt(1.0 + 4.0 * Da / Pe); % auxiliary variable
nom = 2.0 .* exp(0.5 .* Pe .* z) .* (beta .* cosh(0.5 .* Pe .* beta .* (1.0 - z)) + sinh(0.5 .* Pe .* beta .* (1.0 - z)));
denom = (beta^2 + 1.0) * sinh(0.5 * Pe * beta) + 2.0 * beta * cosh(0.5 * Pe * beta);
y = nom ./ denom;
plot(z, y, '--')

% add labels
legend('Numerical', 'Analytical');
xlabel('Dimensionless axial distance [-]')
ylabel('Dimensionless concentration [-]')
grid on;
hold off;