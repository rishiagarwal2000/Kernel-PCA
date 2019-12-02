function [ans] = kernel(xi,xj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = 5;  %%% needs optimization
temp_ans = - dot((xi-xj),(xi-xj));
ans = exp(temp_ans/c^2);
%ans = (dot(xi,xj)+c)^3;
end
