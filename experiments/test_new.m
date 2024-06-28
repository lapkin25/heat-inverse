# See https://github.com/grenkin/joker-fdm

clear all;
#more off;
#format long;

global a
global b
global kappa_a
global alpha
global beta
global gamma
global theta_b1
global theta_b2
global L
global K

# исходные данные
a = 0.92;
b = 18.7;
kappa_a = 0.01;
alpha = 3.333333333;
beta = 10;
gamma = 0.3;
theta_b1 = 0.4;
theta_b2 = 0.4;
L = 50;

# размер расчетной сетки
K = 30;

# источники тепла
global f_fun
global num_sources
f_fun = {@f1, @f2};
num_sources = length(f_fun);

# мощности источников тепла
global q
q = zeros(num_sources);

# известные средние значения температуры по каждому источнику
global r
r = [4, 3];

# источник тепла
function ret = g_fun (x)
  global num_sources
  global f_fun
  global q
  s = 0;
  for i = 1:num_sources
    s += q(i) * f_fun{i}(x);
  endfor
  ret = s;
endfunction

function ret = f_cos(a, b, x)
  if (x >= a && x <= b)
    ret = 0.5*(1 + cos(2 * pi * (x - (a + b) / 2) / (b - a)));
    #if (x < (a + b) / 2)
    #  ret = (x - a) * 2 / (b - a);
    #else
    #  ret = (b - x) * 2 / (b - a);
    #endif
  else
    ret = 0.0;
  endif
endfunction

function ret = f1 (x)
  ret = f_cos(10, 20, x);
#  if (x >= 10 && x <= 20)
#    ret = 1.0;
#  else
#    ret = 0.0;
#  endif
endfunction

function ret = f2 (x)
  ret = f_cos(30, 40, x);
#ret = f_cos(45, 50, x);
#  if (x >= 30 && x <= 40)
#    ret = 1.0;
#  else
#    ret = 0.0;
#  endif
endfunction


function [r_vals, theta] = calc_heat ()
  global a
  global b
  global kappa_a
  global alpha
  global beta
  global gamma
  global theta_b1
  global theta_b2
  global L
  global K
  global f_fun
  global num_sources

  data.N = N = 2;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;
  data.a(2, 1) = alpha;

  fm =  {@(x, p) b*kappa_a*x^4*sign(x),  @(x, p) -b*kappa_a*x;
         @(x, p) -kappa_a*x^4*sign(x), @(x, p) kappa_a*x};
  dfm = {@(x, p) 4*b*kappa_a*abs(x)^3,  @(x, p) -b*kappa_a;
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

  data.b = [ beta, beta ; gamma, gamma ];
  data.w = [ beta * theta_b1, beta * theta_b2 ; gamma * theta_b1 ^ 4, gamma * theta_b2 ^ 4 ];

  data.G = [ ];

  addpath("joker-fdm/bvp1d");
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L, grid_info.nodes);
  data.g = zeros(N, grid_info.nodes);
  data.g(1, :) = arrayfun(@g_fun, xgrid);
  data.g(2, :) = zeros(1, grid_info.nodes);
  for i = 1 : N
    for k = 1 : N
      data.p(i, k, :) = zeros(1, grid_info.nodes);
    endfor
  endfor

  guess = zeros(N, grid_info.nodes);
  tol = 1e-4;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  theta = sol(1, :);
  phi = sol(2, :);

  #{
  figure
  plot(xgrid, theta, "r", "linewidth", 2);
  xlabel("x");
  ylabel("theta");
  figure
  plot(xgrid, phi, "b", "linewidth", 2);
  xlabel("x");
  ylabel("phi");
  #}

  r_vals = zeros(num_sources, 1);
  for i = 1:num_sources
    # вычисляем интеграл от f_fun{i} * theta
    r_vals(i) = trapz(xgrid, arrayfun(f_fun{i}, xgrid) .* theta);
  endfor
endfunction

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

function [A, rhs] = calc_linear_system ()
  global a
  global b
  global kappa_a
  global beta
  global theta_b
  global L
  global K
  global f_fun
  global num_sources

  A = zeros(num_sources);
  rhs = zeros(num_sources, 1);

  data.N = N = 1;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;

  fm =  {@(x, p) 0};
  dfm = {@(x, p) 0};

  data.f = data.df = cell(N, M, N);
  for j = 1 : M
    for i = 1 : N
      for k = 1 : N
        data.f{i, j, k} = fm{i, k};
        data.df{i, j, k} = dfm{i, k};
      endfor
    endfor
  endfor

  data.b = [ beta, beta ];
  data.w = [ 0, 0 ];

  data.G = [ ];

  addpath("joker-fdm/bvp1d");
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L, grid_info.nodes);

  for i = 1 : N
    for k = 1 : N
      data.p(i, k, :) = zeros(1, grid_info.nodes);
    endfor
  endfor

  for eq_ind = 1 : num_sources
    data.g = zeros(N, grid_info.nodes);
    data.g(1, :) = arrayfun(f_fun{eq_ind}, xgrid);

    guess = zeros(N, grid_info.nodes);
    tol = 1e-4;
    sol = solve_bvp1d(grid_info, data, guess, tol);
    eta = sol(1, :);

    s_vals = zeros(num_sources, 1);
    for j = 1:num_sources
      # вычисляем интеграл от f_fun{j} * eta
      s_vals(j) = trapz(xgrid, arrayfun(f_fun{j}, xgrid) .* eta);
    endfor

    A(:, eq_ind) = s_vals;
  endfor

 # data.g = zeros(N, grid_info.nodes);
 # data.w = [ beta * theta_b, beta * theta_b ];

 # guess = zeros(N, grid_info.nodes);
 # tol = 1e-4;
 # sol = solve_bvp1d(grid_info, data, guess, tol);
 # eta = sol(1, :);

 # s_vals = zeros(num_sources, 1);
 # for j = 1:num_sources
 #   # вычисляем интеграл от f_fun{j} * eta
 #   s_vals(j) = trapz(xgrid, arrayfun(f_fun{j}, xgrid) .* eta);
 # endfor

 # rhs = s_vals;
 rhs = zeros(num_sources, 1);
endfunction


#q = [0.1; 0.1; 0.1];
#[r_vals, theta] = calc_heat();

#xgrid = linspace(0, L, 150);
#plot(xgrid, arrayfun(@(x) f1(x) + f2(x) + f3(x), xgrid));

#r_vals
#jac = calc_jacobian_naive(r_vals);
#jac

#jac = calc_jacobian(theta)

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

#[qq, fval, info] = fsolve (@opt_f, q, optimset ("jacobian", "on"))

##q1_min = 0;
##q1_max = 1.5;
##q2_min = 0;
##q2_max = 1.5;
##qn = 10;
##q1grid = linspace(q1_min, q1_max, qn);
##q2grid = linspace(q2_min, q2_max, qn);
##func_vals1 = zeros(qn);
##func_vals2 = zeros(qn);
###j_val = 1;
##q(3) = 0;
##for q1_ind = 1:qn
##  for q2_ind = 1:qn
##    q1_val = q1grid(q1_ind);
##    q2_val = q2grid(q2_ind);
##    q(1) = q1_val;
##    q(2) = q2_val;
##    [r_vals, theta] = calc_heat();
##    func_vals1(q1_ind, q2_ind) = r_vals(1);
##    func_vals2(q1_ind, q2_ind) = r_vals(2);
##  endfor
##endfor
##figure
##contour(q1grid, q2grid, func_vals1', 0:15, 'k', 'ShowText', 'on')
##hold on
##contour(q1grid, q2grid, func_vals2', 0:15, 'k', '--', 'ShowText', 'on')

[B, rhs] = calc_linear_system();


s1_min = -5;
s1_max = 5;
s2_min = -5;
s2_max = 5;
sn = 15;
s1grid = linspace(s1_min, s1_max, sn);
s2grid = linspace(s2_min, s2_max, sn);
func_vals1 = zeros(sn);
func_vals2 = zeros(sn);
#j_val = 1;
q(3) = 0;
for s1_ind = 1:sn
  for s2_ind = 1:sn
    s1_val = s1grid(s1_ind);
    s2_val = s2grid(s2_ind);
    s_val = [s1_val; s2_val];
    q_val = B^-1 * (s_val - rhs);
    q(1) = q_val(1);
    q(2) = q_val(2);
    [r_vals, theta] = calc_heat();
    func_vals1(s1_ind, s2_ind) = r_vals(1);
    func_vals2(s1_ind, s2_ind) = r_vals(2);
  endfor
endfor
figure
contour(s1grid, s2grid, func_vals1', 0:0.2:15, 'k', 'ShowText', 'on')
hold on
contour(s1grid, s2grid, func_vals2', 0:0.2:15, 'k', '--', 'ShowText', 'on')
xlabel("s_1")
ylabel("s_2")

