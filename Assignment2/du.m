function output=du(x,L,k,E)
output=(k*L/(4*pi*pi*E))*(2*pi*k*x*sin(2*pi*k*x/L)+L*cos(2*pi*k*x/L))+(k/(4*pi*pi*E)-4);
end