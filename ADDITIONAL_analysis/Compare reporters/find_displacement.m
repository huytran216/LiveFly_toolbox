function [xout,fun_target] = find_displacement(beta0,fun,y,x,ax,ratio)

    fun_ = subs(fun,y,beta0);
%     fun_inverse = finverse(fun_,x);
%     if numel(fun_inverse)
%         fun_target = arrayfun(@(k) double(subs(fun_,x,k)),ax);
%         xnew = arrayfun(@(k) double(subs(fun_inverse,k)),fun_target/ratio);
%         xout = xnew - ax;
%     else
        % Interpolation:
        open_ax = [ax(1):60];
        open_target = arrayfun(@(k) double(subs(fun_,x,k)),open_ax);
        fun_target = arrayfun(@(k) double(subs(fun_,x,k)),ax);
        xnew = interp1(open_target,open_ax,fun_target/2,'linear','extrap');
        xout = ax - xnew;
%     end