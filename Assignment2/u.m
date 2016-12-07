function output=u(x,L,k,E)
output=(k*L/(4*pi*pi*E))*(L*sin(2*pi*k*x/L)/(k*pi)-x*cos(2*pi*k*x/L))+(k/(4*pi*pi*E)-4)*x+3;
end