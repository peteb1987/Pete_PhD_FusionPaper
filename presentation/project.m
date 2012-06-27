function track = project(t, x2_start, a2)

track = zeros(2, length(t));

track(:,1) = x2_start;

for kk = 2:length(t)
    
    dt = t(kk)-t(kk-1);
    track(1,kk) = track(1,kk-1) + track(2,kk-1)*dt + 0.5*a2*dt^2;
    track(2,kk) = track(2,kk-1) + a2*dt;
    
end

track(2,:) = [];

end