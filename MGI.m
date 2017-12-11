close all;
clear all;

%Définition des données du problème

%l correspond à la longueur du triangle équilatéral représentant le robot

l = 4;

%A correspond aux coordonnées des points situés aux extrémités des jointures.
%On considère que le point A1 se trouve toujours à l'origine du repère.

A1 =[0 0];
A2 = [16 0];
A3 = [8 -5];

A = [A1 ; A2 ; A3];

%G correspond aux coordonnées du centre de gravité du triangle équilatéral
%représentant le robot.

G = [3 ; -1+-4/3];
theta = 0;
RotZ = [cos(theta) -sin(theta);
        sin(theta) cos(theta)];

%-----------------------------------------------------------------------
%----------------------Détermination du MGI-----------------------------

XC1 = [-l/2 ;l/3];
OC1 = G + RotZ * XC1;
l1 = sqrt((A1(1,1)-OC1(1,1))^2 + (A1(1,2) - OC1(2,1))^2);

XC2 = [l/2 ;l/3];
OC2 = G + RotZ * XC2;
l2 = sqrt((A2(1,1)-OC2(1,1))^2 + (A2(1,2) - OC2(2,1))^2);

XC3 = [0 ; -2*l/3];
OC3 = G + RotZ * XC3;
l3 = sqrt((A3(1,1)-OC3(1,1))^2 + (A3(1,2) - OC3(2,1))^2);
