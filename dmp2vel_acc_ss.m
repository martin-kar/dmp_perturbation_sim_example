function [ydot, yddot] = dmp2vel_acc_ss(y0, y, g, tau_adapt, w, x, dt, alpha_e, c, D, alpha_z, beta_z,ya,e,kc,tau)

persistent z
if isempty(z)
    z = 0;
end

psi = exp(-0.5*((x-c).^2).*D);
f  = sum(w'*psi)/sum(psi+1.e-10) * x * (g-y0);
z_dot = 1/tau_adapt * (alpha_z*(beta_z*(g-y) - z) + f);
z = z + dt*z_dot;
ydot = z/tau_adapt;
yddot = (z_dot*tau_adapt-tau*z*2*kc*e*(alpha_e*(ya-y-e)))/tau_adapt^2;

end
