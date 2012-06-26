clup
dbstop if error

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

% Discretise
disc_dt = 2;
disc_t = 1:disc_dt:9;
[~, disc_idx, ~] = intersect(t, disc_t);
disc_x = x(:,disc_idx);
disc_v = v(:,disc_idx);

% Find predicted states
pred_disc_x = disc_x(:,1:end-1)+disc_dt*disc_v(:,1:end-1);
pred_disc_v = disc_v(:,1:end-1);

%% Plot - variable rate model
f1 = figure; hold on
xlim([-1,21]), ylim([-1,8]);
plot(x(1,:), x(2,:), 'linewidth', 2);
plot(disc_x(1,:), disc_x(2,:), 'xr', 'markersize', 10);
h_arrows = arrow(disc_x, disc_x+disc_v, 'length', 10, 'BaseAngle', 30);
set(h_arrows, 'EdgeColor', 'r', 'FaceColor', 'r')
plot(x(1,1), x(2,1), 'g*', 'markersize', 15)
plot(x(1,floor(cp*K)), x(2,floor(cp*K)), 'g*', 'markersize', 15)

export_pdf(f1, 'variable_rate_model.pdf', 16, 12)

%% Plot - HMM
f2 = figure; hold on
xlim([-1,21]), ylim([-1,8]);
% plot(x(1,:), x(2,:), 'linewidth', 2);
plot(disc_x(1,:), disc_x(2,:), 'xr', 'markersize', 10);
h_arrows = arrow(disc_x, disc_x+disc_v, 'length', 10, 'BaseAngle', 30);
set(h_arrows, 'EdgeColor', 'r', 'FaceColor', 'r')
% plot(x(1,1), x(2,1), 'g*', 'markersize', 15)
% plot(x(1,floor(cp*K)), x(2,floor(cp*K)), 'g*', 'markersize', 15)
plot(pred_disc_x(1,:), pred_disc_x(2,:), 'xm', 'markersize', 10);
h_arrows = arrow(pred_disc_x, pred_disc_x+pred_disc_v, 'length', 10, 'BaseAngle', 30);
set(h_arrows, 'EdgeColor', 'm', 'FaceColor', 'm')
for kk = 1:size(pred_disc_x, 2)
    plot_gaussian_ellipsoid(pred_disc_x(:,kk), eye(2));
end

export_pdf(f2, 'HMM_model.pdf', 16, 12)

%% Particle filter demo
% h_vrpf = vrpf_demo;
% 
% for kk = 1:length(h_vrpf)
%     export_pdf(h_vrpf(kk), ['vrpf_demo_' num2str(kk) '.pdf']);
% end

%% Deficiency

f3 = figure; hold on;
xlim([0 2]), ylim([0 2])
dt = 1;
plot_gaussian_ellipsoid([1 1], [dt^3/3 dt^2/2; dt^2/2 dt], 1)
plot_gaussian_ellipsoid([1 1], [dt^3/3 dt^2/2; dt^2/2 dt], 0.5)
xlabel('position'), ylabel('velocity')
export_pdf(f3, 'full_rank.pdf')

f4 = figure; hold on;
xlim([0 2]), ylim([0 2])
dt = 1;
plot_gaussian_ellipsoid([1 1], [dt^3/4 dt^2/2; dt^2/2 dt], 1)
plot_gaussian_ellipsoid([1 1], [dt^3/4 dt^2/2; dt^2/2 dt], 0.5)
xlabel('position'), ylabel('velocity')
export_pdf(f4, 'deficient.pdf')
