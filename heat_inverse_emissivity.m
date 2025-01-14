# See https://github.com/grenkin/joker-fdm

clear all;
#more off;
#format long;

global a
global b
global kappa_a
global alpha
global beta
global gamma1
global gamma2
global theta_b1
global theta_b2
global theta_init
global L
global T
global K
global tnum

# исходные данные
a = 0.92;
b = 18.7;
kappa_a = 0.01;
alpha = 3.333333333;
beta = 10;
theta_b1 = 0.6;  # попробовать другие значения
theta_b2 = 0.6;
theta_init = 0.4;
L = 50;
T = 1;

# коэффициент отражения излучения
gamma1 = 0.3;
gamma2 = 0.5;

# размер расчетной сетки
K = 30;
tnum = 20;


# известные средние значения температуры на каждом конце
global r
r = [4, 4];

function [r_vals, theta, phi] = calc_heat ()
  global a
  global b
  global kappa_a
  global alpha
  global beta
  global gamma1
  global gamma2
  global theta_b1
  global theta_b2
  global theta_init
  global L
  global K
  global T
  global tnum

  tau = T / tnum;

  data.N = N = 2;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;
  data.a(2, 1) = alpha;

  fm =  {@(x, p) b*kappa_a*x^4*sign(x) + x / tau,  @(x, p) -b*kappa_a*x;
         @(x, p) -kappa_a*x^4*sign(x), @(x, p) kappa_a*x};
  dfm = {@(x, p) 4*b*kappa_a*abs(x)^3 + 1 / tau,  @(x, p) -b*kappa_a;
         @(x, p) -kappa_a*4*abs(x)^3, @(x, p) kappa_a};

  data.f = data.df = cell(N, M, N);
  for j = 1 : M
    for i = 1 : N
      for k = 1 : N
        data.f{i, j, k} = fm{i, k};
        data.df{i, j, k} = dfm{i, k};
      endfor
    endfor
  endfor

  data.b = [ beta, beta ; gamma1, gamma2 ];
  data.w = [ beta * theta_b1, beta * theta_b2 ; gamma1 * theta_b1 ^ 4, gamma2 * theta_b2 ^ 4 ];

  data.G = [ ];

  addpath("joker-fdm/bvp1d");
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L, grid_info.nodes);

  theta = zeros(K + 1, tnum + 1);
  phi = zeros(K + 1, tnum + 1);
  theta(:, 1) = theta_init;
  phi(:, 1) = 0;  # phi при t = 0 не определено

  for m = 2:tnum + 1
    data.g = zeros(N, grid_info.nodes);
    data.g(1, :) = theta(:, m - 1) / tau;
    data.g(2, :) = zeros(1, grid_info.nodes);
    for i = 1 : N
      for k = 1 : N
        data.p(i, k, :) = zeros(1, grid_info.nodes);
      endfor
    endfor
    guess = zeros(N, grid_info.nodes);
    guess(1, :) = theta(:, m - 1);
    guess(2, :) = phi(:, m - 1);
    tol = 1e-4;
    sol = solve_bvp1d(grid_info, data, guess, tol);
    theta(:, m) = sol(1, :);
    phi(:, m) = sol(2, :);
  endfor

  r_vals = zeros(1, 2);
  tgrid = linspace(0, T, tnum + 1);
  # вычисляем интеграл от theta
  r_vals(1) = trapz(tgrid, theta(1, :));
  r_vals(2) = trapz(tgrid, theta(K + 1, :));
endfunction

#{
function u = calc_linearized (theta, ii)
  global a
  global b
  global kappa_a
  global alpha
  global beta
  global gamma
  global theta_b
  global L
  global K
  global f_fun
  global num_sources

  data.N = N = 2;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;
  data.a(2, 1) = alpha;

  fm =  {@(x, p) b*kappa_a*4*abs(p)^3*x,  @(x, p) -b*kappa_a*x;
         @(x, p) -kappa_a*4*abs(p)^3*x, @(x, p) kappa_a*x};
  dfm = {@(x, p) b*kappa_a*4*abs(p)^3,  @(x, p) -b*kappa_a;
         @(x, p) -kappa_a*4*abs(p)^3, @(x, p) kappa_a};
  data.f = data.df = cell(N, M, N);
  for j = 1 : M
    for i = 1 : N
      for k = 1 : N
        data.f{i, j, k} = fm{i, k};
        data.df{i, j, k} = dfm{i, k};
      endfor
    endfor
  endfor

  data.b = [ beta, beta ; gamma, gamma ];
  data.w = [ 0, 0 ; 0, 0 ];

  data.G = [ ];

  addpath("joker-fdm/bvp1d");
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L, grid_info.nodes);
  data.g = zeros(N, grid_info.nodes);
  data.g(1, :) = arrayfun(f_fun{ii}, xgrid);
  data.g(2, :) = zeros(1, grid_info.nodes);
  for i = 1 : N
    for k = 1 : N
      data.p(i, k, :) = theta(:);
    endfor
  endfor

  guess = zeros(N, grid_info.nodes);
  tol = 1e-4;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  u = sol(1, :);
  z = sol(2, :);
endfunction

function jac = calc_jacobian_naive (r_vals)
  global q
  global num_sources
  delta_q = 0.0001;  # шаг дифференцирования

  q0 = q;
  for i = 1:num_sources
    q = q0;
    q(i) += delta_q;
    [new_r_vals, theta] = calc_heat();
    jac(:, i) = (new_r_vals - r_vals) / delta_q;
  endfor
endfunction

function jac = calc_jacobian (theta)
  global q
  global num_sources
  global L
  global K
  global f_fun

  xgrid = linspace(0, L, K + 1);
  for i = 1:num_sources
    u_i = calc_linearized(theta, i);
    for j = 1:num_sources
      jac(j, i) = trapz(xgrid, arrayfun(f_fun{j}, xgrid) .* u_i);
    endfor
  endfor
endfunction
#}

[r_vals, theta, phi] = calc_heat();
r_vals

#{
xgrid = linspace(0, L, K + 1);
figure
plot(xgrid, theta(:, tnum + 1), "r", "linewidth", 2);
xlabel("x");
ylabel("theta");
figure
plot(xgrid, phi(:, tnum + 1), "b", "linewidth", 2);
xlabel("x");
ylabel("phi");
#}


#xgrid = linspace(0, L, 150);
#plot(xgrid, arrayfun(@(x) f1(x) + f2(x) + f3(x), xgrid));

#r_vals
#jac = calc_jacobian_naive(r_vals);
#jac

#jac = calc_jacobian(theta)

#{
function [y, jac] = opt_f (x)
  global num_sources
  global q
  global r
  for i = 1:num_sources
    q(i) = x(i);
  endfor
  for i = 1:num_sources
    printf("%f ", q(i));
  endfor
  printf("\n");
  y = zeros (num_sources, 1);
  for j = 1:num_sources
    [r_vals, theta] = calc_heat();
    for i = 1:num_sources
      y(i) = r_vals(i) - r(i);
    endfor
    if (nargout == 2)
      jac = calc_jacobian(theta);
    endif
  endfor
endfunction
#}

#[qq, fval, info] = fsolve (@opt_f, q, optimset ("jacobian", "on"))


gamma_min = 0.01;
gamma_max = 0.5;
gamma_n = 5;
gamma_grid = linspace(gamma_min, gamma_max, gamma_n);
func_vals1 = zeros(gamma_n);
func_vals2 = zeros(gamma_n);
for gamma1_ind = 1:gamma_n
  for gamma2_ind = 1:gamma_n
    gamma1_val = gamma_grid(gamma1_ind);
    gamma2_val = gamma_grid(gamma2_ind);
    gamma1 = gamma1_val;
    gamma2 = gamma2_val;
    [r_vals, theta, phi] = calc_heat();
    func_vals1(gamma1_ind, gamma2_ind) = r_vals(1);
    func_vals2(gamma1_ind, gamma2_ind) = r_vals(2);
  endfor
endfor
figure
contour(gamma_grid, gamma_grid, func_vals1', 'k', 'ShowText', 'on')
hold on
contour(gamma_grid, gamma_grid, func_vals2', 'k', '--', 'ShowText', 'on')

