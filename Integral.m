function valor = Integral(u, k)

%MASCARA = load('mask_boundary_SQUARE.dat');
%radio = 50;
%sin_theta_0 = 0.04273;

%radio_vector = track_radio_mask(radio, sin_theta_0, MASCARA);
%dlmwrite('Radio_vector_2.txt', radio_vector)

radio_vector = fscanf(fopen('Temporal_radial_vector.txt', 'r'), '%f');

integral_function = @(phi) Real_pattern(u, phi, radio_vector) * exp(- 1i * k * phi);
INT = integral(integral_function, -pi, pi, 'ArrayValued', true, 'RelTol', 1e-3, 'AbsTol', 1e-3);
valor = 0.5 * INT / pi