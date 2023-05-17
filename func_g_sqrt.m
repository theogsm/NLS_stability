function g = func_g_sqrt(s, omega, gamma, a1, a3)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
g = omega*realsqrt(s)-func_f_sqrt(s,gamma,a1,a3);
end