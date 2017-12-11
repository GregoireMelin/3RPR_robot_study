clear all;
close all;
clc;


% ----------------------------------- %
% ------  DEFINITION VARIABLES ------ %
% ----------------------------------- %

%%Calcul symbolique
syms l1 l2 c2 rho1 rho2 rho3 l3 c3 d3 phi beta;

%%Donnees numeriques
%c2=1;
%l2=1;
%c3=0;
%d3=1;
%l3=1;
%beta=90*pi/180;
%rho1=4/5;
%rho2=3/2;
%rho3=3/2;

% --------------------------------------------------------- %
% ------  DEFINITION DES DONNEES D'ENTREE DU SYSTEME ------ %
% --------------------------------------------------------- %

%%Grandeurs articulaires et orientation de la plateforme
param_in=[rho1,rho2,rho3,phi];
%%Coordonnees des extremites de la plateforme
T_points=[ 0 0 ; c2 0 ; c3 d3 ];
%%Longueurs des cotes de la plateforme
T_lengths=[l1,l2,beta];

% --------------------------------------- %
% ------  MISE EN EQUATION DU MGD  ------ %
% --------------------------------------- %

%%Mise en place manuelle
A = [ 2*T_lengths(1)*cos(param_in(4))-2*T_points(2,1) , 2*T_lengths(1)*sin(param_in(4));...
      2*T_lengths(2)*cos(param_in(4)+beta)-2*T_points(3,1) , 2*T_lengths(2)*sin(param_in(4)+beta)-2*T_points(3,2)];
B = [ 2*T_points(2,1)*T_lengths(1)*cos(param_in(4))-T_lengths(1)^2-T_points(2,1)^2+param_in(2)^2-param_in(1)^2;...
      2*T_points(3,2)*T_lengths(2)*sin(param_in(4)+beta)+2*T_points(3,1)*T_lengths(2)*cos(param_in(4)+beta)-T_lengths(2)^2-c3^2-T_points(3,2)^2+param_in(3)^2-param_in(1)^2];
%X=A\B;
X=get_MGD_3RPR(param_in,T_points,T_lengths);
simplify(X);



% ------------------------------------------ %
% ------  RECHERCHE DES SINGULARITES  ------ %
% ------------------------------------------ %
%%ATTENTION : Calcul symbolique necessaire ici
%determinant=det(A);
S = solve(det(A)==0,phi);
simplify(S);