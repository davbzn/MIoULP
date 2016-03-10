function [n] = n_BK7(f)
c_const     = 299792458;                                % m/s
% f is frequency vector in PHz
% x is wavelength vector in m: (m/s / (PHz*1e15) = m)
x = c_const./(f*1e15).*1e6;

% refractive coefficient dependence on lambda n(lambda)
n = sqrt(1+1.03961212./(1-0.00600069867./x.^2)+0.231792344./(1-0.0200179144./x.^2)+1.01046945./(1-103.560653./x.^2));
