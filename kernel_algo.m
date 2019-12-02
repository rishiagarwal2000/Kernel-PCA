clc; clear;
%%%%%%  k is my kernel function
%%%%%%  x is my data(2-d) and O is my non linear transformation
%%%%%%  Original Data
x = zeros(1000,1);
y = zeros(1000,1);
load ./concentric_points_data
x = x - sum(x)/1000;
y = y - sum(y)/1000;
figure;
scatter(x,y,2,col);
title('My random points');
xlabel('X');
ylabel('Y');
%%%%%%  k(xi,xj) is my dot product of O(xi).O(xj)
%%%%%%  form a n*n matrix K where n is size of data,
%%%%%%  with Kij = k(xi,xj)
n = 1000;
K = zeros(n);
for i=1:n
    for j=1:n
        K(i,j) = kernel([x(i),y(i)]',[x(j),y(j)]');
    end
end
one_mat = ones(size(K));
K = K - one_mat*K - K*one_mat + one_mat*K*one_mat;
%%%%%%  Find ai's by solving eigenvector equation
[U,S,V] = svd(K);
%%%%%%  Normalize them
alpha = zeros(n);
S = S/n;
for i = 1:n
   alpha(:,i) = U(:,i)/sqrt(S(i,i));
end
%%%%%%  find score(i,k); where i represents data and k represents component
for i = 1:n
    for k = 1:n
        score(i,k) = dot(alpha(:,k), K(:,i));
    end
end
figure;
scatter(score(:,1),score(:,2),2,col);
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Second Principal Component vs First Principal Component');
%%%%%%  find 1/di^2 for each eiegen vector
beta = 1;  %%%needs optimization
d_2 = zeros(n);
for i = 1:n
    for k = 1:n
        d_2(i,k) = (K(i,i)-2*score(i,k)*beta+beta*beta)^(-1); 
    end
end
%%%%%%  find sum_d_2_k for each component

sum_d_2 = zeros(n,1);
for i = 1:n
    sum_d_2(i) = sum(d_2(:,i));
end
%%%%%%  find gamma(j) for each test point
gamma = zeros(n,1);
for j = 1:n
    for i = 1:n
        gamma(j) = gamma(j) + score(j,i)*alpha(j,i);
    end
end
%%%%%%  find 1/di^2 for most dominant eiegen vector
  %%%needs optimization
d_22 = zeros(n);
for i = 1:n
    for k = 1:n/2
        beta = k;
        d_22(i,k) = (K(i,i)-2*score(i,1)*beta+beta*beta)^(-1);
        d_22(i,k+n/2) = (K(i,i)-2*score(i,2)*beta+beta*beta)^(-1);
    end
end
%%%%%%  find sum_d_2_k for each beta

sum_d_22 = zeros(n,1);
for i = 1:n
    sum_d_22(i) = sum(d_22(:,i));
end
%%%%%%  find denoised points

%%%%%%  Actually these methods do something else, which I need to understand
%%%%%% Method 1(Failed)
%new_x = zeros(n,1);
%new_y = zeros(n,1);
%for i = 1:n
%    new_x(i) = dot(d_2(:,i),x)/sum_d_2(i);
%    new_y(i) = dot(d_2(:,i),y)/sum_d_2(i);
%end
%figure;
%plot(new_x,new_y);
%%%%%% Method 2
%new_x = zeros(n,1);
%new_y = zeros(n,1);
%for i = 1:n
%   for counter = 1:10
%    tempX = new_x(i); tempY = new_y(i);
%    new_x(i) = 0; new_y(i) = 0; total = 0;
%    for j = 1:n
%        t = kernel([tempX,tempY]',[x(j),y(j)]');
%        new_x(i) = new_x(i) + gamma(j)*t*x(j);
%        new_y(i) = new_y(i) + gamma(j)*t*y(j);
%        total = total + gamma(j)*t;
%    end
%    new_x(i) = new_x(i)/total;
%    new_y(i) = new_y(i)/total;
%   end
%end
%figure;
%scatter(new_x,new_y);
%%%%%%% Method 3
%new_x = zeros(n,1);
%ew_y = zeros(n,1);
%for i = 1:n
%    new_x(i) = dot(d_22(:,i),x)/sum_d_22(i);
%    new_y(i) = dot(d_22(:,i),y)/sum_d_22(i);
%end
%figure;
%scatter(new_x,new_y);
%%%%%%% Method 4
%new_x = U(:,1)'*K';
%new_y = U(:,2)'*K';
%figure;
%scatter(new_x,new_y);
