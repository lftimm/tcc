// CONFIGURAÇÕES
Djet = 1.0;
Zjet = 4*Djet;
splitHeight = 2;
dz = Djet/10;

// GEOMETRIA
// Face da frente
Point(1) = {0,0,0,0};
Point(2) = {0,Zjet,0,0};
Point(3) = {0.5*Djet,Zjet,0,0};
Point(4) = {0.5*Djet,15*Djet,0,0};
Point(5) = {20*Djet,15*Djet,0,0};
Point(6) = {20*Djet,0,0,0};
Point(7) = {20*Djet,Zjet,0,0};
Point(8) = {0.5*Djet,0,0,0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 6};
Line(7) = {6, 8};
Line(8) = {8, 1};
Line(9) = {8, 3};
Line(10) = {3, 7};

Curve Loop(1) = {1, 2, -9, 8};
Plane Surface(1) = {1};

Curve Loop(2) = {10, 6, 7, 9};
Plane Surface(2) = {2};

Curve Loop(3) = {5, -10, 3, 4};
Plane Surface(3) = {3};

Transfinite Surface {1} = {1, 2, 3, 8};
Transfinite Surface {2} = {8, 3, 7, 6};
Transfinite Surface {3} = {7, 5, 4, 3};


prog=1.005;
num_eles=250;
Transfinite Curve {1, 9} = num_eles+1 Using Progression prog;
Transfinite Curve {6} = num_eles+1 Using Progression 1/prog;

Transfinite Curve {10, 4} = 750+1 Using Progression 1.0020;
Transfinite Curve {7} = 750+1 Using Progression 1/1.0020;

Transfinite Curve {8, 2} = 50+1 Using Progression 1;

prog2=1.005;
Transfinite Curve {3} = 130+1 Using Progression prog2;
Transfinite Curve {5} = 130+1 Using Progression 1/prog2;

Recombine Surface {1,2,3};

Extrude {0, 0, 0.1} {
  Surface{1}; Surface{2}; Surface{3}; Layers {1}; Recombine;
}

Physical Surface("Radius", 4) = {19};
Physical Surface("Impinging Jet", 1) = {23};
Physical Surface("Wall", 2) = {31, 49};
Physical Surface("Slip Wall", 3) = {71};
Physical Surface("Top Pressure Outlet", 8) = {75};

Physical Surface("Pressure Outlet", 6) = {45,63};
Physical Volume("Fluid", 7) = {1, 2, 3, 4, 5};

Physical Surface("Extrusion", 5) = {3,76,54,2,71,19,63,45,23,15,32,49,75,31,1};


