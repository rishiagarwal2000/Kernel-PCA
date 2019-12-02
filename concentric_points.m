clc;clear;
x = zeros(1000,1);
y = zeros(1000,1);

for i=1:1000
    r1 = rand(1);
    if r1<0.33
        r = 0.5*rand(1);
    elseif r1<0.66
        r = 2+rand(1);
    else
        r = 4+rand(1);
    end
    theta = 2*pi*rand(1);
    x(i) = r*cos(theta);
    y(i) = r*sin(theta);
end
