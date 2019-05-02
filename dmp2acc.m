function [ yddot] = dmp2acc(y0, y, ydot, g, tau, w, x, c, D, alpha_z, beta_z)

psi = exp(-0.5*((x-c).^2).*D);
f  = sum(w'*psi)/sum(psi+1.e-10) * x * (g-y0);
yddot = 1/tau^2 * (alpha_z*(beta_z*(g-y) - tau*ydot) + f);

end

