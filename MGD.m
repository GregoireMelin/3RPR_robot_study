clear all;
close all;
clc;
% ------------------------------------------------ %
% ---- Instanciation des variables --------------- %
% ------------------------------------------------ %
syms x y t;



%Coordonnees des points du triangle
c2 = 15.91;
c3 = 0;
d3 = 10;
T_points=[0 0;c2 0;c3 d3];

%Parametres articulaires (mis au carre)
rho1 = 14.98^2;
rho2 = 15.38^2;
rho3 = 12^2;
p_joint=[rho1 rho2 rho3];
theta = 0.882603;

%Longueurs des cotes de la plateforme
l2 = 17.04;
l3 = 20.84;
l1 = l2^2+l3^2 - 2 * l2 * l3 * cos(theta); %Al kashi

T_lengths=[l1 l2 l3];

%Necessaire pour un changement de variable pour la resolution
sinus_phi = (2*t)/(1+t^2);
cosinus_phi= (1-t^2)/(1+t^2);

% ------------------------------------------------ %
% ---- MISE EN PLACE DU MGD ---------------------- %
% ------------------------------------------------ %

%On exprime rho1, rho2 et rho3 en fonction des longueurs des plateformes
%rho1^2=x^2 + y^2;  [1]
%rho2^2=(x+l2*cos(phi)-c2)^2 + (y+l2*sin(phi))^2;   [2]
%rho^3=(x+l3*cos(phi+beta)-c3)^2+(y+l3*sin(phi+beta)-d3)^2; [3]
%En faisant [3] - [1] et [2] - [1] : on obtient le systeme suivant

R = 2*l2*cosinus_phi-2*T_points(2,1);
S = 2*l2*sinus_phi;
Q = -2*T_points(2,1)*l2*cosinus_phi + l2^2 + T_points(2,1)^2;
U = 2*l3*(cosinus_phi*cos(theta)-sinus_phi*sin(theta))-2*T_points(3,1);
V = 2*l3*(sinus_phi*cos(theta)+cosinus_phi*sin(theta)) -2*T_points(3,2);
W = -2*T_points(3,2)*l3*(sinus_phi*cos(theta)+cosinus_phi*sin(theta)) - 2*T_points(3,1)*l3*(cosinus_phi*cos(theta)-sinus_phi*sin(theta)) + l3^2 + T_points(3,1)^2 + T_points(3,2)^2;

eqn_xy = (S*(p_joint(3) - p_joint(1) - W) - V*(p_joint(2) - p_joint(1) - Q))^2 +(R*(p_joint(3) - p_joint(1) - W) - U*(p_joint(2) - p_joint(1) - Q))^2 -p_joint(1)*(R*V - S*U)^2;
eqn_xy = simplify(eqn_xy);

% ------------------------------------------------ %
% ------------- RESOLUTION  ---------------------- %
% ------------------------------------------------ %

%On cherche a annuler le polynome d'ou le fait de n'etudier que son denominateur
P = numden(eqn_xy);

%Mise sous forme de polynome
C = coeffs(P(1,1));
C=fliplr(C); % Necessaire

%Calcul des racines du polynome
rotz = roots(C);

%On a p_positione phi = 2*atan(t)*180/pi : Calcul des orientations p_positionsibles pour
%des coordonnees articulaires donnees
rotz_phi = atan(rotz)*2*180/3.14;

for i=1:6
    p_position(i,1) = -(S*(p_joint(3) - p_joint(1) - W) - V*(p_joint(2) - p_joint(1) - Q))/(R*V - S*U);
    p_position(i,2) =  (R*(p_joint(3) - p_joint(1) - W) - U*(p_joint(2) - p_joint(1) - Q))/(R*V - S*U);
    p_position(i,1) = subs( p_position(i,1) ,t, rotz(i,1));
    p_position(i,2) = subs( p_position(i,2) ,t, rotz(i,1));
    p_position(i,3) = rotz(i,1);
    difference=p_joint(1)-(p_position(i,1)^2+p_position(i,2)^2); % Doit etre nulle si les p_positionitions sont justes
end
