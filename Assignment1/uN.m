function output=uN(x,he,N)
i=floor(x*N)+1;
output=a(i)*theta1(x,he,N)+a(i+1)*theta2(x,he,N);
end