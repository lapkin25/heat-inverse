# Попытка реализовать метод "предиктор-корректор" в нестационарном варианте
# Данная реализация работает странно

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
global tnum
global delta_t
global q_init
global theta_init

# исходные данные
a = 0.92;
b = 18.7;
kappa_a = 0.01;
alpha = 3.333333333;
#beta = 10;
gamma = 0.3;
#gamma = beta / a * alpha
beta = gamma / alpha * a
theta_b1 = 0.4;
theta_b2 = 0.4;
theta_init = 0.4;
q_init = 0;
L = 50;

global lambda
lambda = 1.0;

# размер расчетной сетки
K = 50;
tnum = 300;
delta_t = 1.0;

# источники тепла
global f_fun
global num_sources
f_fun = {@f1, @f2, @f3};
num_sources = length(f_fun);

# мощности источников тепла
global q
q = zeros(num_sources, 1);

# заданные средние значения температуры по каждому источнику
global r
r = [4, 4, 3];

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
    #ret = 0.5*(1 + cos(2 * pi * (x - (a + b) / 2) / (b - a)));
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

function ret = f3 (x)
  ret = f_cos(45, 50, x);
#  if (x >= 45 && x <= 50)
#    ret = 1.0;
#  else
#    ret = 0.0;
#  endif
endfunction


function A = calc_linear_system ()
  global a
  global b
  global kappa_a
  global beta
  global L
  global K
  global f_fun
  global num_sources

  A = zeros(num_sources);

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
endfunction


function phi = calc_phi (theta)
  global alpha
  global b
  global kappa_a
  global gamma
  global theta_b1
  global theta_b2
  global L
  global K
  global f_fun
  global num_sources

  A = zeros(num_sources);

  data.N = N = 1;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = alpha;

  fm =  {@(x, p) kappa_a*x};
  dfm = {@(x, p) kappa_a};

  data.f = data.df = cell(N, M, N);
  for j = 1 : M
    for i = 1 : N
      for k = 1 : N
        data.f{i, j, k} = fm{i, k};
        data.df{i, j, k} = dfm{i, k};
      endfor
    endfor
  endfor

  data.b = [ gamma, gamma ];
  data.w = [ gamma * theta_b1 ^ 4, gamma * theta_b2 ^ 4 ];

  data.G = [ ];

  addpath("joker-fdm/bvp1d");
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L, grid_info.nodes);

  for i = 1 : N
    for k = 1 : N
      data.p(i, k, :) = zeros(1, grid_info.nodes);
    endfor
  endfor

  data.g = zeros(N, grid_info.nodes);
  data.g(1, :) = kappa_a * theta .^ 4;

  guess = zeros(N, grid_info.nodes);
  tol = 1e-4;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  phi = sol(1, :);
endfunction


function [q_vals, r_vals] = calc_heat ()
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
  global q
  global r
  global tnum
  global delta_t
  global q_init
  global theta_init
  global lambda

  phi_init = calc_phi(theta_init * ones(K + 1, 1));
  Mat = calc_linear_system();

  r_vals = zeros(num_sources, tnum + 1);
  rad_vals = zeros(num_sources, tnum + 1);
  q_vals = zeros(num_sources, tnum + 1);

  data.N = N = 2;  # число уравнений
  data.M = M = 1;  # число слоев
  data.a(1, 1) = a;
  data.a(2, 1) = alpha;

  fm =  {@(x, p) b*kappa_a*x^4*sign(x) + x / delta_t,  @(x, p) -b*kappa_a*x;
         @(x, p) -kappa_a*x^4*sign(x), @(x, p) kappa_a*x};
  dfm = {@(x, p) 4*b*kappa_a*abs(x)^3 + 1 / delta_t,  @(x, p) -b*kappa_a;
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

  for j = 1:num_sources
    q_vals(j, 1) = q_init;
    r_vals(j, 1) = trapz(xgrid, arrayfun(f_fun{j}, xgrid) .* theta_init);
    rad_vals(j, 1) = trapz(xgrid, arrayfun(f_fun{j}, xgrid) .* phi_init);
    q(j) = q_init;
  endfor

  theta_cur = theta_init * ones(1, grid_info.nodes);

  for tt = 2 : tnum + 1
    data.g = zeros(N, grid_info.nodes);
    data.g(1, :) = arrayfun(@g_fun, xgrid) + 1 / delta_t * theta_cur;
    data.g(2, :) = zeros(1, grid_info.nodes);
    for i = 1 : N
      for k = 1 : N
        data.p(i, k, :) = zeros(1, grid_info.nodes);
      endfor
    endfor

    guess = zeros(N, grid_info.nodes);
    guess(1, :) = theta_cur;
    tol = 1e-4;
    sol = solve_bvp1d(grid_info, data, guess, tol);
    theta = sol(1, :);
    phi = sol(2, :);

    theta_cur = theta;

    for i = 1:num_sources
      # вычисляем интеграл от f_fun{i} * theta
      r_vals(i, tt) = trapz(xgrid, arrayfun(f_fun{i}, xgrid) .* theta);
      rad_vals(i, tt) = trapz(xgrid, arrayfun(f_fun{i}, xgrid) .* phi);
    endfor

    #for i = 1:num_sources
    q_vals(:, tt) = q_vals(:, tt - 1) + lambda * delta_t * Mat^-1 * (r' - r_vals(:, tt)); #(r' - r_vals(:, tt - 1));
    # здесь в теории должно быть +=, а не -=
    q_vals(:, tt) -= b * alpha / a * Mat^-1 * (rad_vals(:, tt) - rad_vals(:, tt - 1));
    q = q_vals(:, tt);
    #endfor

  endfor

endfunction


[q_vals, r_vals] = calc_heat()

t_grid = (0:tnum) * delta_t;
figure
plot(t_grid, r_vals(1, :), 'color', 'black', "-", t_grid, r_vals(2, :), ...
  'color', 'black', "-.", t_grid, r_vals(3, :), 'color', 'black', "--")
legend("1", "2", "3")
xlabel("t")
ylabel("Average temperature")

figure
plot(t_grid, q_vals(1, :), 'color', 'black', "-", t_grid, q_vals(2, :), ...
  'color', 'black', "-.", t_grid, q_vals(3, :), 'color', 'black', "--")
legend("1", "2", "3")
xlabel("t")
ylabel("Intensities of heat sources")

