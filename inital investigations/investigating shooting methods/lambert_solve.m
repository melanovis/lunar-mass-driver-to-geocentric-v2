function [v_1,v_2,delta_t] = lambert_solve(r1, r2, dt_normalised, mu)

coord = norm(r2-r1);
semiparameter = (norm(r2)+norm(r1)+coord)/2;

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
u_1 = r1/norm(r1);
u_2 = r2/norm(r2);
A = sqrt(mu/(4*a_n))*cot(alpha_solve/2);
B = sqrt(mu/(4*a_n))*cot(beta_solve/2);
u_c = (r2-r1)/coord;
v_1 = (B+A)*u_c + (B-A)*u_1;
v_2 = (B+A)*u_c - (B-A)*u_2;


end