function [param_out,A,B] = get_MGD_3RPR(param_in,T_points,T_lengths,phi)
   
   %param_in=[rho1,rho2,rho3,phi];
   %T_points=[ 0 0 ; c2 0 ; c3 d3 ];
   %T_lengths=[l1,l2,beta];

   A = [ 2*T_lengths(1)*cos(phi)-2*T_points(2,1) , 2*T_lengths(1)*sin(phi);...
         2*T_lengths(2)*cos(phi+T_lengths(3))-2*T_points(3,1) , 2*T_lengths(2)*sin(phi+T_lengths(3))-2*T_points(3,2)];

   B = [ 2*T_points(2,1)*T_lengths(1)*cos(phi)-T_lengths(1)^2-T_points(2,1)^2+param_in(2)^2-param_in(1)^2;...
         2*T_points(3,2)*T_lengths(2)*sin(phi+T_lengths(3))+2*T_points(3,1)*T_lengths(2)*cos(phi+T_lengths(3))-T_lengths(2)^2-T_points(3,1)^2-T_points(3,2)^2+param_in(3)^2-param_in(1)^2];

   param_out=A\B;

end