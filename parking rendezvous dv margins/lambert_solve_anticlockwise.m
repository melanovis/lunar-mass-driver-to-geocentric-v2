function [v_1,v_2,delta_t] = lambert_solve_anticlockwise(r_1, r_2, dt_normalised, mu)

coord = norm(r_2-r_1);
semiparameter = (norm(r_2)+norm(r_1)+coord)/2;

a_min = semiparameter/2;
a_max = 2*semiparameter;

alpha_solve = 2*asin(sqrt(semiparameter/(2*a_min)));
beta_solve = 2*asin(sqrt((semiparameter-coord)/(2*a_min)));

tmin = (sqrt(2)/3)*sqrt((semiparameter^3)/mu)*(1-((semiparameter-coord)/semiparameter)^(3/2));
tmax = sqrt((a_min^3)/mu)*(alpha_solve-beta_solve-(sin(alpha_solve)-sin(beta_solve)));

t_target = interp1([0,1],[tmin,tmax], dt_normalised);
t_prev = tmax;

if semiparameter == coord
    semiparameter = semiparameter + 1e-3;
end

for n=1:50
    a_n = (a_min+a_max)/2;
    alpha_solve = 2*asin(sqrt(semiparameter/(2*a_n)));
    beta_solve = 2*asin(sqrt((semiparameter-coord)/(2*a_n)));

    delta_t = sqrt((a_n^3)/mu)*(alpha_solve-beta_solve-(sin(alpha_solve)-sin(beta_solve)));
    
    if delta_t < t_target
        a_max = a_n;
    else
        a_min = a_n;
    end

    if abs(delta_t-t_prev) < 1e-7
        break
    end
    t_prev = delta_t;
end

%finding v 
u_1 = r_1/norm(r_1);
u_2 = r_2/norm(r_2);
A = sqrt(mu/(4*a_n))*cot(alpha_solve/2);
B = sqrt(mu/(4*a_n))*cot(beta_solve/2);
u_c = (r_2-r_1)/coord;
v_1 = (B+A)*u_c + (B-A)*u_1;
v_2 = (B+A)*u_c - (B-A)*u_2;


%we're only considering anticlockwise v_1 vectors

if dot(cross(r_1,v_1),[0,0,1]) < 0
    v_1 = -v_1;
    v_2 = -v_2;
    transfer_period_full = sqrt(((4*pi^2)/mu)*a_n^3);
    delta_t = transfer_period_full - delta_t;
end


end