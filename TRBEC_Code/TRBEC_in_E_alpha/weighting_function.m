%computes the weight function as a function of alpha_eq, used for
%interpolation so no inputs required, output is in radians!
function [W,alpha_eq_mid] = weighting_function()

dalpha_eq = 1/100;%DO NOT TAKE TOO SMALL A STEP IN alpha_eq (causes issues when taking derivative of Y(y))
dtheta = 1/100000;%DO NOT TAKE A BIG STEP IN dtheta (causes issues when taking derivative of Y(y))

alpha_eq = (deg2rad(1):dalpha_eq:deg2rad(179))';

alpha_eq_length = length(alpha_eq);
y = sin(alpha_eq);

%need to go through each alpha_eq one at time, calculate theta_m
%and integrate over theta in order to obtain a y value,
%then once all y's have been obtained take derivative and are left with

%will need to compute T(y) for each alpha_eq bin, and for each time
%step since the loss cone changes over time
I_3 = zeros(alpha_eq_length,1);
for i = 1:alpha_eq_length
    theta_m = colatitude_m(alpha_eq(i));
    theta = theta_m:dtheta:pi/2;
    
    Y_y_theta = sin(theta).*sqrt(1+3.*(cos(theta)).^2).*((1-(y(i).^2.*sin(theta).^(-6).*sqrt(1+3*cos(theta).^2))).^(1/2));
    
    %take integral
    I_3(i) = (trapz(theta,Y_y_theta,2));
    
    %if Nan value is given (will occur at 90 deg) then replace with 0 (the actual value)
    if isnan(I_3(i))
       I_3(i) = 0; 
    end
end

%take derivative (taken from matlab documentation where derivative is approximated as diff()/dh)
I_3_derivative = diff(I_3)./dalpha_eq;
%derivative is estiamted at the mid points, not the original bin boundaries
alpha_eq_mid = (alpha_eq(2:end)+alpha_eq(1:end-1))/2;
%interpolate I_3 at the mid points (so it's the same size as the
%derivative)
I_3 = interp1(alpha_eq,I_3,alpha_eq_mid);

%calculate weighting function
W = abs(real(cos(alpha_eq_mid).*sin(alpha_eq_mid).*I_3 - sin(alpha_eq_mid).^2.*I_3_derivative));

%explicitly define W(0) = W(pi) = 0
alpha_eq_mid = [0;alpha_eq_mid;pi];
W = [0;W;0];
end