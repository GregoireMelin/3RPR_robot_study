clear all;
close all;
clc;

% ----------------------------------- %
% ------  DEFINITION VARIABLES ------ %
% ----------------------------------- %

%%Calcul symbolique
syms rho1 rho2 rho3;
%%Donnees numeriques
c2=1;
l2=1;
l1=1;
c3=0;
d3=1;
l3=1;
beta=60*pi/180;

rho1=1.4142;
rho2=11.0454;
rho3=5;

% --------------------------------------------------------- %
% ------  DEFINITION DES DONNEES D'ENTREE DU SYSTEME ------ %
% --------------------------------------------------------- %

%%Grandeurs articulaires et orientation de la plateforme
param_in=[rho1,rho2,rho3];
%%Coordonnees des extremites de la plateforme
T_points=[ 0 0 ; c2 0 ; c3 d3 ];
%%Longueurs des cotes de la plateforme
T_lengths=[l1,l2,beta];
% -------------------------------------------------------- %
% ------  CALCUL DE L'ORIENTATION DE LA PLATEFORME  ------ %
% -------------------------------------------------------- %
poly_phi=[T_points(3,1)*(param_in(1)^2-param_in(2)^2+4*T_points(2,1)^2-4*T_points(2,1)*T_points(3,1))+T_points(2,1)*(param_in(3)^2-param_in(1)^2), ...
          T_points(3,2)*(8*T_points(2,1)*T_points(3,1)-4*T_points(2,1)^2+param_in(2)^2-param_in(1)^2), ...
          T_points(3,1)*(param_in(1)^2-param_in(2)^2)+param_in(3)^2*T_points(2,1)-4*T_points(3,2)^2*T_points(2,1)-param_in(1)^2*T_points(2,1), ...
          T_points(3,2)*(param_in(2)^2-param_in(1)^2)];

sol_poly_phi=roots(poly_phi);

% Dans le polynome x = tan(phi/2) d'ou phi = 2* atan2(x)
sol_poly_phi(1)=2*atan(sol_poly_phi(1));
sol_poly_phi(2)=2*atan(sol_poly_phi(2));
sol_poly_phi(3)=2*atan(sol_poly_phi(3));
%phi=get_phi(param_in,T_points);

%3 solutions pour phi source de singularite
phi1=sol_poly_phi(1);
phi2=sol_poly_phi(2);
phi3=sol_poly_phi(3);
phi=phi1;
% --------------------------------------- %
% ------  MISE EN EQUATION DU MGD  ------ %
% --------------------------------------- %

R=2*T_lengths(1)*cos(phi)-2*T_points(2,1);
S=2*T_lengths(1)*sin(phi);
Q= 2*T_points(2,1)*T_lengths(1)*cos(phi)-T_lengths(1)^2-T_points(2,1)^2+param_in(2)^2-param_in(1)^2;

U=2*T_lengths(2)*cos(phi+beta)-2*T_points(3,1);
V=2*T_lengths(2)*sin(phi+beta)-2*T_points(3,2);
W= 2*T_points(3,2)*T_lengths(2)*sin(phi+T_lengths(3))+2*T_points(3,1)*T_lengths(2)*cos(phi+T_lengths(3))-T_lengths(2)^2-T_points(3,1)^2-T_points(3,2)^2+param_in(3)^2-param_in(1)^2;

A = [R , S ; ...
     U , V];
B = [Q ; W];
%Cas ou RS - SU != 0
X=A\B;
determinant=det(A);




