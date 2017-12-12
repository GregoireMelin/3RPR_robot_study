function [poly_phi,sol_poly_phi] = get_phi(param_in,T_points)
   
   %param_in=[rho1,rho2,rho3,phi];
   %T_points=[ 0 0 ; c2 0 ; c3 d3 ];
   %T_lengths=[l1,l2,beta];
    poly_phi=[T_points(3,1)*(param_in(1)^2-param_in(2)^2+4*T_points(2,1)^2-4*T_points(2,1)*T_points(3,1))+T_points(2,1)*(param_in(3)^2-param_in(1)^2), ...
              T_points(3,2)*(8*T_points(2,1)*T_points(3,1)-4*T_points(2,1)^2+param_in(2)^2-param_in(1)^2), ...
          T_points(3,1)*(param_in(1)^2-param_in(2)^2)+param_in(3)^2*T_points(2,1)-4*T_points(3,2)^2*T_points(2,1)-param_in(1)^2*T_points(2,1), ...
          T_points(3,2)*(param_in(2)^2-param_in(1)^2)];

    sol_poly_phi=roots(poly_phi);
    % Dans le polynome x = tan(phi/2) d'ou phi = 2* atan2(x)
    sol_poly_phi(1)=2*atan(sol_poly_phi(1));
    sol_poly_phi(2)=2*atan(sol_poly_phi(2));
    sol_poly_phi(3)=2*atan(sol_poly_phi(3));

end