function output=duN(x,J,a,N)
he=2*J;
i=floor(x*N)+1;
theta=(2*x-(2*i-1)*he)/he;
output=(a(2*i+1)*(theta+0.5)+a(2*i)*(-2*theta)+a(2*i-1)*(theta-0.5))/J;
end