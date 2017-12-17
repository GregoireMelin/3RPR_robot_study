%-----------------------------------------------------------------------
%----------------------Détermination du MGD-----------------------------

%Paramètres d'entrée

%Les articulations du robot sont reliés à chaque point du triangle A1A2A3, où A1 se situe à
%l'origine, et A2 a pour ordonnée zéro. Chaque point de ce triangle est
%relié à un point du triangle R1R2R3, triangle représentant le robot.

A1 =[0 0];
A2 = [15.91 0];
A3 = [0 10];

c2 = A2(1,1);
c3 = A3(1,1);
d3 = A3(1,2);

l2 = 17.04;        % longueur A1A2
l3 = 20.84;        %longueur A1A3
theta = 0.882603;  %angle en radian entre A1A2 et A1A3


p1 = 14.98^2;      %longueur de l'articulation au carré reliant A1 à R1
p2 = 15.38^2;      %longueur de l'articulation au carré reliant A2 à R2
p3 = 12^2;         %longueur de l'articulation au carré reliant A3 à R3

G = zeros(6,3); %La matrice G va stocker la position du centre de gravité du triangle R1R2R3 (position en abscisse au sein de la première colonne, position en ordonnée dans la seconde), ainsi que la rotation en z du triangle (au sein de la troisième colonne) en degrée

%--------------------------------------------------------------------------
%------------------------Calcul du MGD-------------------------------------

syms t;                 % variable représentant la tangente de la moitié de l'angle de rotation en z du centre de gravité de R1R2R3 (point G)

sphi = (2*t)/(1+t^2);   % variable représentant le sinus de l'angle de rotation en z (angle phi)
cphi= (1-t^2)/(1+t^2);  % variable représentant le cosinus de l'angle de rotation en z (angle phi)

R = 2*l2*cphi-2*c2;
S = 2*l2*sphi;
Q = -2*c2*l2*cphi + l2^2 + c2^2;

U = 2*l3*(cphi*cos(theta)-sphi*sin(theta))-2*c3;
V = 2*l3*(sphi*cos(theta)+cphi*sin(theta)) -2*d3;
W = -2*d3*l3*(sphi*cos(theta)+cphi*sin(theta)) - 2*c3*l3*(cphi*cos(theta)-sphi*sin(theta)) + l3^2 + c3^2 + d3^2;


% R1x = -(SA1 - VA2)/(RV - SU)
% R1y = (RA2 - UA2)/(RV - SU)
% R1x et R1y sont exprimés en fonction de t


A1 = p3 - p1 - W;
A2 = p2 - p1 - Q;

finalEq = (S*A1 - V*A2)^2 +(R*A1 - U*A2)^2 -p1*(R*V - S*U)^2;    %les valeurs de t annulant finalEq correspondent aux valeurs de t pouvant être prise par le système avec les données d'entrées fournies

finalEq = simplify(finalEq);

P = numden(finalEq);     %le polynôme caractéristique du MGD correspond au numérateur obtenu. Il est extrait du quotient par la fonction numden

C = coeffs(P(1,1));      %extraction des coefficient du polynome caractéristique

C = fliplr(C);           %inversion de la matrice (les coefficient du polynome sont inscrit dans le sens croissant de leur monôme correspondant par la fonction coeffs, fliplr permet de placer les coefficient selon un ordre décroissant)

tValues = roots(C);         %détermination des racines du polynome. Ces racines correspond à la valeur de l'angle t pouvant être prise, en radian

rotz = atan(tValues)*2;      %détermination de l'angle rotz pouvant être pris, en radian

G(1:6,3) = rotz*180/3.14;  %stockage de la rotation du triangle sur l'axe z en degrée

for i=1:6                          %détermination des coordonnées de R1 R2 R3

     R1(i,1) = -(S*A1 - V*A2)/(R*V - S*U);
     R1(i,2) =  (R*A1 - U*A2)/(R*V - S*U);

     R1(i,1) = subs(R1(i,1) ,t, tValues(i,1));
     R1(i,2) = subs(R1(i,2) ,t, tValues(i,1));


     R2(i,1) = R1(i,1) +l2*cos(rotz(i,1));
     R2(i,2) = R1(i,2) +l2*sin(rotz(i,1));

     R3(i,1) = R1(i,1) + l3*cos(rotz(i,1)+theta);
     R3(i,2) = R1(i,2) + l3*sin(rotz(i,1)+theta);

     %détermination du centre de gravité de R1R2R3
     G(i,1) = (R1(i,1) + R2(i,1) + R3(i,1))/3;
     G(i,2) = (R1(i,2) + R2(i,2) + R3(i,2))/3;

end

%G est la matrice stockant les différentes solutions du MGD
%Chaque ligne de cette matrice donne une solution au système
%La première colonne indique la coordonnée en abscisse de chaque point
%La seconde colonne indique la coordonnée en ordonnée de chaque point
%La troisième colonne indique la rotation sur l'axe z du triangle
