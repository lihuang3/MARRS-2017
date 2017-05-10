N = 100;
x = zeros(N-1,1);
y = x;
A = 10000;
p = 0.01;
n = 500;
c2 = 1e-5;
c1 = 1e-4;
for ii =2:N
    r = 1/ii;
    x(ii-1) = ii;
    %y(ii-1) = c1*A*p*(1-a)/(b-2)*((b-1)^(1+log(A)/log(b))-1);
    y(ii-1) = c1*p*A*(1-(1/r-1)^(log(A)./log(1./r)+1))./(2-1./r);
    y1(ii-1) = c2*n*A*((1-r)^2-(1-r)^(2*log(A)/log(1/r)+1))/(1-(1-r)^2);

    
end

figure(1)
plot(x,y);hold on
plot(x,y1,'-r'); hold off
