function value = Complex_pattern(u, Phi, radio_vector)
    n = 6; % Ceros de Bessel.

    zeros = [1.219669891266504491 2.233130594381528944 3.238315484166236047 4.241062863796069671 5.242764376870182375 6.243921689864488478 7.244759868719957474 8.245394913952042160 9.245892684949467366 10.24629334875491615 11.24662279487788474];

    % Raíces complexas.

    u_complex = [0.5967 1.7837 3.6420 4.3039 5.2119];
    v_complex = [0.5225 0.5268 0 0 0];
    
    while (abs(Phi) > pi)
        if (Phi > pi)
            Phi = Phi - 2 * pi;
        else
            Phi = Phi + 2 * pi;
        end
    end
    
    m = ceil(1800 + round(10 * Phi * 180 / pi));
    u = u * radio_vector(m+1);
    h_1 = 1; h_2 = 1; h_f = 1;

    for j = 1:(n-1)
        h_1 = h_1 * u_complex(j)^4 + v_complex(j)^4 + u^4 + 2 * u_complex(j)^2 * v_complex(j)^2 + 2 * u^2 * (v_complex(j)^2 - u_complex(j)^2);
        h_2 = h_2 * u_complex(j)^2 + v_complex(j)^2;
        h_f = h_f * (1 - (u / zeros(j))^2);
    end
    
    value = (2 * besselj(1, pi * u) / (h_f * pi * u))^2 * h_1 / h_2^2;
end