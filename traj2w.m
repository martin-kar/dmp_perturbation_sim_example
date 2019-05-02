function w = traj2w(y_demo, dt, tau, c, D, alpha_z, beta_z, alpha_x, n_kernel)

P = length(y_demo); % Nbr time steps
g = y_demo(end);
y0 = y_demo(1);
ydot_demo = [diff(y_demo); 0]/dt;
yddot_demo = [diff(ydot_demo); 0]/dt;

f_target = tau^2*yddot_demo -alpha_z*(beta_z*(g-y_demo) -tau*ydot_demo);


x = zeros(size(y_demo));
x(1) = 1;
for t=2:P
    x_dot = -alpha_x*x(t-1)/tau;
    x(t) = x(t-1) + x_dot*dt;
end

psi = zeros(P, n_kernel);


for t = 1:P
    exp(-0.5*((x(t)-c).^2).*D);
    psi(t,:) = exp(-0.5*((x(t)-c).^2).*D);
    psi(t,:);
end

s = x*(g-y0);

w = zeros(n_kernel,1);
for i = 1:n_kernel
    Gamma_i = diag(psi(:,i));
    w(i) = s'*Gamma_i * f_target / (s'*Gamma_i*s);  
end
