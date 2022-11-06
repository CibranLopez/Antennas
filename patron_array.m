function result = patron_array(u,v)
    global Imn_plot Imn_xy_exc_x Imn_xy_exc_y
    F = 0;
    % CASO SIN SUBARRAY
    %     spacing_sub = 0.5;
    %     M = 5;
    %     N = 5;
    rho = sqrt(u.^2+v.^2);
    if rho <= 1
        for n = 1:numel(Imn_plot)
            if abs(Imn_plot(n)) >= 2e-2
                %abs(Imn_plot(n))
                arg_exp=2*pi.*(Imn_xy_exc_x(n).*u+Imn_xy_exc_y(n).*v);
                F = F + Imn_plot(n).*exp(1i.*arg_exp); 
            end
        end
        % SIN SUBARRAY
        %F = F .* f_subarray(M,N,spacing_sub,u,v);
    end
    result = abs(F).^2;
end