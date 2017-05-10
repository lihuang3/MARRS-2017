N = 50;
x = zeros(N,1);
y = x;
y1 = x;
A = 1e4;
p = 0.01;
n = 100;
c2 = 1e-5;
c1 = 1e-3;
for ii =1:N
    r = 1-1/(N+1)*ii;
    x(ii) = r;
    %y(ii-1) = c1*A*p*(1-a)/(b-2)*((b-1)^(1+log(A)/log(b))-1);
    y1(ii) = c2*n*A*((1-r)^2-(1-r)^(2*log(A)/log(1/r)+1))/(1-(1-r)^2);

    
end

figure(1);clf

plot(x,y1,'-r');

figure(2);clf

M = 50;
x = zeros(M,1);
y = x;

for ii = 1:M
    x(ii) = ii;
   
    for jj = ii:-1:1
        r = (jj-1)/jj;
        y(ii) = y(ii) + A*n/ii*c2*(1-r)^1%;(1*(ii-jj+1)); 
    end
end
plot(x,y)