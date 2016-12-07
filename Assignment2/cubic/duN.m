function output=duN(x,J,a,N)
he=2*J;
i=floor(x*N)+1;
theta=(2*x-(2*i-1)*he)/he;
output=(a(3*i+1)*((9/16)*(3*theta*theta+2*theta-1/9))+a(3*i)*((-27/16)*(3*theta*theta+theta*2/3-1))+a(3*i-1)*((27/16)*(3*theta*theta-theta*2/3-1))+a(3*i-2)*((-9/16)*(3*theta*theta-2*theta-1/9)))/J;
end