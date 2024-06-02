
%given alpha_eq returns theta_m, input and output are in radians!
function theta_interp = colatitude_m(alpha_interp)
%since theta_m is symetric around 90, colat(alpha) = colat(180-alpha)
if alpha_interp > pi/2
    alpha_interp = pi - alpha_interp;
end
theta_m = 0:0.01:pi/2;
alpha_eq = asin(sqrt(((sin(theta_m).^6))./(sqrt(3*(cos(theta_m)).^2 + 1))));
theta_interp = interp1(alpha_eq,theta_m,alpha_interp);
end