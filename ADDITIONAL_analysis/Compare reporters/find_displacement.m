function [xout,fun_target] = find_displacement(beta0,fun,y,x,ax,ratio)

    fun = subs(fun,y,beta0);
    fun_inverse = finverse(fun,x);
    
    fun_target = arrayfun(@(k) double(subs(fun,x,k)),ax);
    xnew = arrayfun(@(k) double(subs(fun_inverse,k)),fun_target/ratio);
    xout = xnew - ax;
    
    
    