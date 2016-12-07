function output=uN(x,he,a,N)
i=floor(x*N)+1;
theta=(2*x-(2*i-1)*he)/he;
output=a(2*i-1)*thetahat1(theta)+a(2*i)*thetahat2(theta)+a(2*i+1)*thetahat3(theta);
end