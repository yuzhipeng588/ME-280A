function output=uN(x,he,a,N)
i=floor(x*N)+1;
theta=((2*x-(2*i-1)*he)/he);
output=a(3*i-2)*thetahat1(theta)+a(3*i-1)*thetahat2(theta)+a(3*i)*thetahat3(theta)+a(3*i+1)*thetahat4(theta);
end