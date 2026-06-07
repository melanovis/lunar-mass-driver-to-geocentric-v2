function state_out = threebody_ODE(t, state_in, earth_mass, lunar_mass)

    G = 6.6743e-11;

    r_rocket = state_in(1:3).';
    v_rocket = state_in(4:6).';
    r_earth = state_in(7:9).';
    v_earth = state_in(10:12).';
    r_luna = state_in(13:15).';
    v_luna = state_in(16:18).';

    r_earthluna = r_earth - r_luna;
    r_unit_earthluna = norm(r_earthluna);

    acc_earth =  -G * lunar_mass * r_earthluna / r_unit_earthluna^3;
    acc_luna = G * earth_mass * r_earthluna / r_unit_earthluna^3;
    
    acc_rocket = (-G*earth_mass)*(r_rocket-r_earth)/(norm(r_rocket-r_earth)^3) + (-G*lunar_mass)*(r_rocket-r_luna)/(norm(r_rocket-r_luna)^3);

    state_out = [v_rocket, acc_rocket, v_earth, acc_earth, v_luna, acc_luna].';
end