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

% Proposse some initial particles
pts_a2 = normrnd(0, 1, [5,1]);

pts_v2 = 

