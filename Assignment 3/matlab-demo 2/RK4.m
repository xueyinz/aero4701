function xnext = RK4( f ,x  ,dt) 
    
    k1 = f(x);
    k2 = f(x+k1*dt/2);
    k3 = f(x+k2*dt/2);
    k4 = f(x+k3*dt); 
    
    xnext = x + dt*(k1+2*k2+2*k3+k4)/6; 
    
end