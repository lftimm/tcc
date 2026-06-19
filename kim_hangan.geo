//DEFINIÇÕES
dz = 0.1;
Djet = 1;
Zjet = 4*Djet;
y = 0.5*Zjet;

slipWallYElements = 30;
slipWallXElements = 30;

wallYElements = 50;
wallXElements = 100;

jetYElements = 25;

//PONTOS
Point(1) = {0,0,0,1.0};
Point(2) = {0,Zjet,0,1.0};
Point(3) = {0.5*Djet,Zjet,0,1.0};
Point(4) = {0.5*Djet,10*Djet,0,1.0};
Point(5) = {10*Djet,10*Djet,0,1.0};
Point(6) = {10*Djet,0,0,1.0};
Point(7) = {0,y,0,0,1.0};
Point(9) = {10*Djet,y,0,1.0};
Point(10) = {10*Djet,y,0,1.0};
Point(11) = {10*Djet,Zjet,0,1.0};

//GEOMETRIA
Line(1) = {1, 7};
Line(2) = {7, 9};
Line(3) = {9, 6};
Line(4) = {6, 1};
Line(5) = {7, 2};
Line(6) = {2, 11};
Line(7) = {11, 9};
Line(8) = {11, 3};
Line(9) = {3, 4};
Line(10) = {4, 5};
Line(11) = {5, 11};

//SUPERFICIES
Curve Loop(1) = {9, 10, 11, 8};
Plane Surface(1) = {1};

Curve Loop(2) = {6, 7, -2, 5};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 4, 1, 2};
Plane Surface(3) = {3};

//DISCRETIZAÇÃO
Transfinite Surface {1} = {3, 4, 5, 11};
Transfinite Surface {2} = {7, 2, 11, 9};
Transfinite Surface {3} = {1, 7, 9, 6};

Transfinite Curve {1, 3} = wallYElements Using Progression 1;
Transfinite Curve {5, 7} = jetYElements Using Progression 1;
Transfinite Curve {4, 2, 6} = wallXElements Using Progression 1;

Transfinite Curve {9, 11} = slipWallYElements+1 Using Progression 1;
Transfinite Curve {8, 10} = slipWallXElements+1 Using Progression 1;

Recombine Surface {3, 2, 1};

//EXTRUSÃO
Extrude {0, 0, dz} {
  Surface{1}; Surface{2}; Surface{3}; 
}
