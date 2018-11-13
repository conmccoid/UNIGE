function y = rk4(f,tspan,y)

h = mean(tspan(2:end) - tspan(1:end-1));

for t = tspan(1:end-1)
    
    k1 = f(t,y);
    k2 = f(t + h/2, y + h*k1/2);
    k3 = f(t + h/2, y + h*k2/2);
    k4 = f(t + h, y + h*k3);
    
    y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

end