close all;
clear all;

%  coordonnees globales ---->[ MGI ]-----> coordonnees articulaires

%%Donnees d'entree : coordonnees globales
%Position

%Points du robot lies au sol
A =[0 0];
C = [16 0];
F = [8 -5];
points_ground=[A,C,F];

%longueurs du triangle (cas du triangle equilateral l1=l2=l3)
l= 4;

%Soit G le centre de gravite du triangle
G = [3 ;(-1+-4/3)];

%Orientation
%Soit phi, l'orientation de la plateforme
phi = 0;
rot_z = [cos(phi) -sin(phi) ; sin(phi) cos(phi)]; %Matrice permettant le changement de repere

%% Calcul des coordonnees articulaires resultantes
% On calcule les coordonnees des points dans le repere de la plateforme mobile ayant pour origine le
%centre de gravite (R1) et on calcule ensuite les memes coordonnees dans le repere lie au sol ayant pour
%origine le point A (R0).

%Calcul des coordonnees des points dans le repere R0
B_R1 = [-l/2 ; l/3];                                     % coordonnees de B dans R1
B = G + rot_z * B_R1;                                   % coordonnees de B dans R0
D_R1 = [l/2 ; l/3];                                     % coordonnees de D dans R1
D = G + rot_z * D_R1;                                   % coordonnees de D dans R0
E_R1 = [0 ; -2*l/3];                                     % coordonnees de E dans R1
E = G + rot_z * E_R1;                                    % % coordonnees de E dans R0

%Calcul des coordonnes articulaires
rho1 = sqrt((A(1,1)-B(1,1))^2 + (A(1,2) - B(2,1))^2); % norme du vecteur AB = rho1
rho2 = sqrt((C(1,1)-D(1,1))^2 + (C(1,2) - D(2,1))^2); % norme du vecteur CD = rho2
rho3 = sqrt((F(1,1)-E(1,1))^2 + (F(1,2) - E(2,1))^2);  % norme du vecteur EF = rho3
joint_coordinates = [rho1,rho2,rho3]; %Longueur des articulations consultable par consultation de la variable joint_coordinates.
