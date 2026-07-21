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
#beta = 10;  # 0.05;
#beta = 0.05;
gamma = 0.3;
beta = gamma / alpha * a; #0.05;  #0.08; #0.08; #10;
theta_b1 = 0.4;
theta_b2 = 0.8;
L = 50;

# размер расчетной сетки
K = 70;

# источники тепла
global f_fun
global h_fun
global num_sources
f_fun = {@f1, @f2};
h_fun = {@h1, @h2};
num_sources = length(f_fun);

# мощности источников тепла
global q
q = zeros(num_sources);

# известные средние значения температуры по каждому источнику
global r
r = [1.8, 1.8];
#r = [4.8, 2.1];

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
  #  ret = 0.5*(1 + cos(2 * pi * (x - (a + b) / 2) / (b - a)));
    if (x < (a + b) / 2)
      ret = (x - a) * 2 / (b - a);
    else
      ret = (b - x) * 2 / (b - a);
    endif
  else
    ret = 0.0;
  endif
endfunction

function ret = f1 (x)
  ret = f_cos(10, 20, x);
endfunction

function ret = f2 (x)
  ret = f_cos(35, 45, x);
endfunction

function ret = h1 (x)
  ret = f_cos(20, 35, x); #- 2
endfunction

function ret = h2 (x)
  ret = f_cos(35, 40, x); #- 2
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
  global h_fun
  global num_sources

  data.N = N = 2;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;
  data.a(2, 1) = alpha;

  fm =  {@(x, p) b*kappa_a*x^4*sign(x),  @(x, p) -b*kappa_a*x;
         @(x, p) -kappa_a*x^4*sign(x), @(x, p) kappa_a*x};
  dfm = {@(x, p) 4*b*kappa_a*abs(x)^3,  @(x, p) -b*kappa_a;
         @(x, p) -kappa_a*4*abs(x)^3, @(x, p) kappa_a};
#  fm =  {@(x, p) b*kappa_a*abs(x)^5*sign(x),  @(x, p) -b*kappa_a*x;
#         @(x, p) -kappa_a*abs(x)^5*sign(x), @(x, p) kappa_a*x};
#  dfm = {@(x, p) 5*b*kappa_a*abs(x)^4,  @(x, p) -b*kappa_a;
#         @(x, p) -kappa_a*5*abs(x)^4, @(x, p) kappa_a};

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
  tol = 1e-6;
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
    # вычисляем интеграл от h_fun{i} * theta
    r_vals(i) = trapz(xgrid, arrayfun(h_fun{i}, xgrid) .* theta);
  endfor
endfunction

function [A, rhs] = calc_linear_system ()
  global a
  global b
  global kappa_a
  global beta
  global theta_b1
  global theta_b2
  global L
  global K
  global f_fun
  global h_fun
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

  data.b = [ [beta, beta] ];
  data.w = [ [0, 0] ];

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
    tol = 1e-6;
    sol = solve_bvp1d(grid_info, data, guess, tol);
    eta = sol(1, :);

    s_vals = zeros(num_sources, 1);
    for j = 1:num_sources
      # вычисляем интеграл от f_fun{j} * eta
      s_vals(j) = trapz(xgrid, arrayfun(h_fun{j}, xgrid) .* eta);
    endfor

    A(:, eq_ind) = s_vals;
  endfor

  data.g = zeros(N, grid_info.nodes);
  data.w = zeros(1, 2);
  data.w = [[ beta * theta_b1, beta * theta_b2 ]];

  guess = zeros(N, grid_info.nodes);
  tol = 1e-6;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  eta = sol(1, :);

  s_vals = zeros(num_sources, 1);
  for j = 1:num_sources
    # вычисляем интеграл от h_fun{j} * eta
    s_vals(j) = trapz(xgrid, arrayfun(h_fun{j}, xgrid) .* eta);
  endfor

  rhs = s_vals;
endfunction

function [A, rhs] = calc_linear_system_2 ()
  global a
  global alpha
  global b
  global kappa_a
  global beta
  global gamma
  global theta_b1
  global theta_b2
  global L
  global K
  global f_fun
  global h_fun
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

  bc_coef = min(beta / a, gamma / alpha) * a;

  #data.b = [ [beta, beta] ];
  data.b = [ [bc_coef, bc_coef] ];
  data.w = [ [0, 0] ];

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
    tol = 1e-6;
    sol = solve_bvp1d(grid_info, data, guess, tol);
    eta = sol(1, :);

    s_vals = zeros(num_sources, 1);
    for j = 1:num_sources
      # вычисляем интеграл от h_fun{j} * eta
      s_vals(j) = trapz(xgrid, arrayfun(h_fun{j}, xgrid) .* eta);
    endfor

    A(:, eq_ind) = s_vals;
  endfor

  data.g = zeros(N, grid_info.nodes);
  data.w = zeros(1, 2);
  #data.w = [[ beta * theta_b1, beta * theta_b2 ]];
  data.w = [[ bc_coef * theta_b1, bc_coef * theta_b2 ]];


  guess = zeros(N, grid_info.nodes);
  tol = 1e-6;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  eta = sol(1, :);

  s_vals = zeros(num_sources, 1);
  for j = 1:num_sources
    # вычисляем интеграл от h_fun{j} * eta
    s_vals(j) = trapz(xgrid, arrayfun(h_fun{j}, xgrid) .* eta);
  endfor

  rhs = s_vals;
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
  tol = 1e-6;
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
  global h_fun

  xgrid = linspace(0, L, K + 1);
  for i = 1:num_sources
    u_i = calc_linearized(theta, i);
    for j = 1:num_sources
      jac(j, i) = trapz(xgrid, arrayfun(h_fun{j}, xgrid) .* u_i);
    endfor
  endfor
endfunction

q = [0; 0];
[r_vals, theta] = calc_heat();

[B, rhs] = calc_linear_system_2();  # система координат с уменьшенным beta

# начальное приближение для оптимальной системы координат
s_init = [B(1, 1), B(1, 2), B(2, 1), B(2, 2)];

function p = p_c(y)
  c = 100000;
  if (y <= 0)
    p = y ^ 2;
  else
    p = c * y ^ 2;
  endif
endfunction

function y = opt_f_slow (s_coef)
  global q

  S = zeros(2);
  S(1, 1) = s_coef(1);
  S(1, 2) = s_coef(2);
  S(2, 1) = s_coef(3);
  S(2, 2) = s_coef(4);

  sum = 0.0;
  qn = 3;
  q_max = 3;
  q1grid = linspace(0, q_max, qn);
  q2grid = linspace(0, q_max, qn);
  for q1_ind = 1:qn
    for q2_ind = 1:qn
      q1_val = q1grid(q1_ind);
      q2_val = q2grid(q2_ind);
      q = [q1_val; q2_val]
      [r_vals, theta] = calc_heat();
      M = calc_jacobian(theta);
      # TODO: предрассчитать якобиан во всех узлах коллокации, чтобы не повторять расчет
      # TODO: сгенерировать узлы коллокации с помощью Latin Hypercube Sampling
      grad_F1 = M(1, :)';
      grad_F2 = M(2, :)';
      dF1_dsk = S^-1 * grad_F1;
      dF2_dsk = S^-1 * grad_F2;
      dF1_ds2 = dF1_dsk(2)
      dF2_ds1 = dF2_dsk(1)
      sum += p_c(dF1_ds2) + p_c(dF2_ds1);
    endfor
  endfor

  y = sum * det(S);

  S
  y
endfunction


function y = opt_f(s_coef)
  global collocation_nodes_num
  global grad_F1_collocation
  global grad_F2_collocation

  S = zeros(2);
  S(1, 1) = s_coef(1);
  S(1, 2) = s_coef(2);
  S(2, 1) = s_coef(3);
  S(2, 2) = s_coef(4);

  sum = 0.0;
  for i = 1:collocation_nodes_num
    grad_F1 = grad_F1_collocation(i, :)';
    grad_F2 = grad_F2_collocation(i, :)';
    dF1_dsk = S^-1 * grad_F1;
    dF2_dsk = S^-1 * grad_F2;
    dF1_ds2 = dF1_dsk(2);
    dF2_ds1 = dF2_dsk(1);
    sum += p_c(dF1_ds2) + p_c(dF2_ds1);
  endfor

  y = sum * det(S);

  S
  y
endfunction


# узлы коллокации
# в каждом столбце матрицы collocation_nodes координаты узла
qn = 3;
global collocation_nodes_num
collocation_nodes_num = qn * qn;
collocation_nodes = zeros(collocation_nodes_num, 2);

q_max = 0.5;
q1grid = linspace(0, q_max, qn);
q2grid = linspace(0, q_max, qn);
cnt = 0;
for q1_ind = 1:qn
  for q2_ind = 1:qn
    q1_val = q1grid(q1_ind);
    q2_val = q2grid(q2_ind);
    cnt += 1;
    collocation_nodes(cnt, :) = [q1_val, q2_val];
  endfor
endfor

# расчет градиента в узлах коллокации...
global grad_F1_collocation
global grad_F2_collocation
grad_F1_collocation = zeros(collocation_nodes_num, 2);
grad_F2_collocation = zeros(collocation_nodes_num, 2);

for i = 1:collocation_nodes_num
  i
  q1_val = collocation_nodes(i, 1);
  q2_val = collocation_nodes(i, 2);
  q = [q1_val; q2_val]
  [r_vals, theta] = calc_heat();
  Mat = calc_jacobian(theta);
  grad_F1_collocation(i, :) = Mat(1, :);
  grad_F2_collocation(i, :) = Mat(2, :);
endfor


# раскомментировать для запуска оптимизации

[s_coef, fval, info] = fsolve(@opt_f, s_init)
opt_S = [s_coef(1), s_coef(2); s_coef(3), s_coef(4)]


# Построение графика...

s1_min = 150;
s1_max = 1000;
s2_min = 150;
s2_max = 1000;
sn = 5;
s1grid = linspace(s1_min, s1_max, sn);
s2grid = linspace(s2_min, s2_max, sn);
func_vals1 = zeros(sn);
func_vals2 = zeros(sn);

for s1_ind = 1:sn
  for s2_ind = 1:sn
    s1_val = s1grid(s1_ind);
    s2_val = s2grid(s2_ind);
    s_val = [s1_val; s2_val];
    #q_val = opt_S^-1 * (s_val - rhs);
    q_val = opt_S^-1 * s_val;
    q(1) = q_val(1);
    q(2) = q_val(2);

    if (q(1) >= 0 && q(2) >= 0)
      q
      [r_vals, theta] = calc_heat();
      func_vals1(s1_ind, s2_ind) = r_vals(1);
      func_vals2(s1_ind, s2_ind) = r_vals(2);
      theta_positive = true;
      for i = 1:length(theta)
        if (theta(i) < 1e-6)
          theta_positive = false;
        endif
      endfor

      #Mat = calc_jacobian(theta);
      #angle1(s1_ind, s2_ind) = atan2(Mat(1, 2), Mat(1, 1));
      #angle2(s1_ind, s2_ind) = atan2(Mat(2, 2), Mat(2, 1));
    else
      func_vals1(s1_ind, s2_ind) = NaN;
      func_vals2(s1_ind, s2_ind) = NaN;
    endif
  endfor
endfor

figure
contour(s1grid, s2grid, func_vals1', 0:0.05:15, 'k', 'ShowText', 'on')
hold on
contour(s1grid, s2grid, func_vals2', 0:0.05:15, 'k', '--', 'ShowText', 'on')
xlabel("s_1")
ylabel("s_2")



# вывести якобиан и значения производных dF1/ds2, dF2/ds1 во всех узлах коллокации

sum = 0.0;
for i = 1:collocation_nodes_num
  q1_val = collocation_nodes(i, 1);
  q2_val = collocation_nodes(i, 2);
  q = [q1_val; q2_val]
  grad_F1 = grad_F1_collocation(i, :)';
  grad_F2 = grad_F2_collocation(i, :)';
  dF1_dsk = opt_S^-1 * grad_F1;
  dF2_dsk = opt_S^-1 * grad_F2;
  dF1_ds1 = dF1_dsk(1)
  dF1_ds2 = dF1_dsk(2)
  dF2_ds1 = dF2_dsk(1)
  dF2_ds2 = dF2_dsk(2)
  sum += p_c(dF1_ds2) + p_c(dF2_ds1);
endfor
y = sum * det(opt_S)

