clup;
dbstop if error

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

%% Generate a trajectory
dt = 0.01;
t = 0:dt:10;
K = length(t);

a = zeros(2,K);
v = zeros(2,K);
x = zeros(2,K);

v(:,1) = [2, 2]';
x(:,1) = [0, 0]';

cp = 0.6;
a(2,1:floor(cp*K)) = -0.5;
a(2,floor(cp*K)+1:end) = 1;

for kk = 2:K
    v(:,kk) = v(:,kk-1) + dt*a(:,kk-1);
    x(:,kk) = x(:,kk-1) + dt*v(:,kk-1) + 0.5*dt^2*a(:,kk-1);
end

% Discretise it
disc_dt = 2;
disc_t = 1:disc_dt:9;
[~, disc_idx, ~] = intersect(t, disc_t);
disc_x = x(:,disc_idx);
disc_v = v(:,disc_idx);

%% Ok. Now treat disc_t, disc_x and disc_v as observations, and mock-pf it
t = disc_t;
y = [disc_x; disc_v];

% Propose some initial particles
Np = 10;
pts_a1 = normrnd(0, 0.1, [Np,1]);
pts_a2 = normrnd(-0.5, 0.5, [Np,1]);

%% t=1

% Project them forward
tt = 0:dt:1; Kt = length(tt);
pts_x2 = zeros(Np,Kt);
for ii = 1:Np
    pts_x2(ii,:) = project(tt, [0 2]', pts_a2(ii));
end

% x1 coords
pts_x1 = zeros(Np,Kt);
for ii = 1:Np
    pts_x1(ii,:) = project(tt, [0 2]', pts_a1(ii));
end
% pts_x1 = repmat(project(tt, [0 2]', 0), Np, 1);

% Plot it
figure, hold on,
xlim([-1,21]), ylim([-1,8]);
plot(y(1,:), y(2,:), 'rx', 'markersize', 10);
plot(pts_x1', pts_x2', 'b')

%% t=2

% Project them forward
tt = 0:dt:3; Kt = length(tt);
pts_x2 = zeros(Np,Kt);
for ii = 1:Np
    pts_x2(ii,:) = project(tt, [0 2]', pts_a2(ii));
end

% x1 coords
pts_x1 = zeros(Np,Kt);
for ii = 1:Np
    pts_x1(ii,:) = project(tt, [0 2]', pts_a1(ii));
end
% pts_x1 = repmat(project(tt, [0 2]', 0), Np, 1);

% Plot it
figure, hold on,
xlim([-1,21]), ylim([-1,8]);
plot(y(1,:), y(2,:), 'rx', 'markersize', 10);
plot(pts_x1', pts_x2', 'b')

best = [5,6]; rest = [1:4, 7:10];

% Plot resampled
figure, hold on,
xlim([-1,21]), ylim([-1,8]);
plot(y(1,:), y(2,:), 'rx', 'markersize', 10);
plot(pts_x1(best,:)', pts_x2(best,:)', 'b')

% RM
pts_a2(rest) = normrnd(-0.5, 0.1, [8,1]);
for ii = rest
    pts_x2(ii,:) = project(tt, [0 2]', pts_a2(ii));
end

% Plot again
figure, hold on,
xlim([-1,21]), ylim([-1,8]);
plot(y(1,:), y(2,:), 'rx', 'markersize', 10);
plot(pts_x1', pts_x2', 'b')

%% t=3

% Project them forward
tt = 0:dt:5; Kt = length(tt);
pts_x2 = zeros(Np,Kt);
for ii = 1:Np
    pts_x2(ii,:) = project(tt, [0 2]', pts_a2(ii));
end

% x1 coords
pts_x1 = zeros(Np,Kt);
for ii = 1:Np
    pts_x1(ii,:) = project(tt, [0 2]', pts_a1(ii));
end
% pts_x1 = repmat(project(tt, [0 2]', 0), Np, 1);

% Plot it
figure, hold on,
xlim([-1,21]), ylim([-1,8]);
plot(y(1,:), y(2,:), 'rx', 'markersize', 10);
plot(pts_x1', pts_x2', 'b')