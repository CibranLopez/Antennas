function value = Real_pattern(u, phi_rad, radio_vector)
    n = 6; % Ceros de Bessel.
    M = 2; % Número de ceros en v.
    zeros = [1.219669891266504491 2.233130594381528944 3.238315484166236047 4.241062863796069671 5.242764376870182375 6.243921689864488478 7.244759868719957474 8.245394913952042160 9.245892684949467366 10.24629334875491615 11.24662279487788474];

    % Raíces complexas.
    
    u_real = [1.0225 3.0445 1.0225 3.0445 5.4319 6.1570 7.1293];
    v_real = [1.1424 1.0829 -1.1424 -1.0829 0.0000 0.0000 0.0000];
    
    phi_deg = phi_rad * 180 / pi;
    m = ceil(1+((180+round(10*phi_deg)/10)/0.1));
    u = u * radio_vector(m);    
    num = 1; den = 1;
    
    for k=1:n+M-1
        num = num.*(1-u^2/(u_real(k)+1i.*v_real(k)).^2);
    end
    
    for k=1:n+M-1        
        den = den.*(1-u^2./zeros(k)^2);
    end
    
    value = 2*besselj(1,pi.*u)/(pi.*u).*num./den;
end