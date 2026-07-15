// CONFIGURAÇÕES
dz = 0.1;
Djet = 0.0381;
Zjet = 4*Djet;
splitHeight = 0.45*Zjet;

// Aqui usei este valor como o governante no tamanho da malha
// Defini mais ou menos o tamanho da região onde ocorrem os vórtices
// Usei proporções entre os diferentes elementos pra compor uma malha mais uniforme
elementCoef = 350;

// GEOMETRIA
// Face da frente
Point(1) = {0,0,0,1.0};
Point(2) = {0,Zjet,0,1.0};
Point(3) = {0.5*Djet,Zjet,0,1.0};
Point(4) = {0.5*Djet,10*Djet,0,1.0};
Point(5) = {10*Djet,10*Djet,0,1.0};
Point(6) = {10*Djet,0,0,1.0};
Point(7) = {0,splitHeight,0,0,1.0};
Point(9) = {10*Djet,splitHeight,0,1.0};
Point(10) = {10*Djet,splitHeight,0,1.0};
Point(11) = {10*Djet,Zjet,0,1.0};
Point(12) = {0.5*Djet,splitHeight,0,1.0};
Point(13) = {0.5*Djet,0,0,1.0};

Line(1) = {1, 7};
Line(2) = {7, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 11};
Line(7) = {11, 9};
Line(8) = {9, 6};
Line(11) = {3, 11};
Line(13) = {12, 13};
Line(14) = {3, 12};
Line(15) = {12, 7};
Line(16) = {9, 12};
Line(17) = {6, 13};
Line(18) = {13, 1};

Curve Loop(1) = {2, 3, 14, 15};
Plane Surface(1) = {1};

Curve Loop(2) = {1, -15, 13, 18};
Plane Surface(2) = {2};

Curve Loop(3) = {16, 13, -17, -8};
Plane Surface(3) = {3};

Curve Loop(4) = {7, 16, -14, 11};
Plane Surface(4) = {4};

Curve Loop(5) = {4, 5, 6, -11};
Plane Surface(5) = {5};

Transfinite Surface {5} = {3, 4, 5, 11};
Transfinite Surface {4} = {12, 3, 11, 9};
Transfinite Surface {3} = {13, 12, 9, 6};
Transfinite Surface {2} = {1, 7, 12, 13};
Transfinite Surface {1} = {7, 2, 3, 12};

// Retirada de dados da malha
p13[] = Point{13};
p12[] = Point{12};
p6[] = Point{6};
p1[] = Point{1};

line13Length = p12[1] - p13[1];
line17Length = p6[0]-p13[0];
line18Length = p13[0]-p1[0];
line1718Proportion = line17Length/line18Length;
plane3SideProportions = line17Length/line13Length;


elementCoefElementsX = elementCoef;
elementCoefElementssplitHeight = elementCoefElementsX/plane3SideProportions;

// Dividindo as linhas, aqui estipulei 10% pra suavizar a transição entre as regiões
Transfinite Curve {1, 13, 8} = elementCoefElementssplitHeight+1 Using Progression 1;
Transfinite Curve {2, 7, 14} = elementCoefElementssplitHeight*1.10+1 Using Progression 1;
Transfinite Curve {6, 4} = elementCoefElementssplitHeight*1.1*1.1+1  Using Progression 1;

Transfinite Curve {18, 15, 3} = elementCoefElementsX/line1718Proportion+1 Using Progression 1;
Transfinite Curve {17, 16, 11, 5} = elementCoefElementsX+1 Using Progression 1;

// Recombinando
Recombine Surface {5};
Recombine Surface {4};
Recombine Surface {3};
Recombine Surface {2};
Recombine Surface {1};

// Extrusão
Extrude {0, 0, dz} {
  Surface{5}; Surface{4}; Surface{3}; Surface{2}; Surface{1}; Layers {1}; Recombine;
}

// Grupos Físicos
// Tolerância para a captura geométrica
tol = 0.01;

// Identificação das faces geradas pela extrusão nas fronteiras (Baseadas nas coordenadas X, Y, Z)
surf_wall[] = Surface In BoundingBox{-tol, -tol, -tol, 10*Djet+tol, tol, dz+tol};
surf_slip[] = Surface In BoundingBox{0.5*Djet-tol, Zjet-tol, -tol, 0.5*Djet+tol, 10*Djet+tol, dz+tol};
surf_outlet_vert[] = Surface In BoundingBox{10*Djet-tol, -tol, -tol, 10*Djet+tol, 10*Djet+tol, dz+tol};
surf_outlet_hori[] = Surface In BoundingBox{0.5*Djet-tol, 10*Djet-tol, -tol, 10*Djet+tol, 10*Djet+tol, dz+tol};

surf_jet[] = Surface In BoundingBox{-tol, Zjet-tol, -tol, 0.5*Djet+tol, Zjet+tol, dz+tol};
surf_axis[] = Surface In BoundingBox{-tol, -tol, -tol, tol, Zjet+tol, dz+tol};

// Identificação dos planos transversais (Z = 0 e Z = dz) para forçar o escoamento 2D
planos_z0[] = Surface In BoundingBox{-tol, -tol, -tol, 10*Djet+tol, 10*Djet+tol, tol};
planos_zdz[] = Surface In BoundingBox{-tol, -tol, dz-tol, 10*Djet+tol, 10*Djet+tol, dz+tol};
todos_volumes[] = Volume "*";

Physical Surface("Wall", 129) = {surf_wall[]};
Physical Surface("Slip Wall", 130) = {surf_slip[]};
Physical Surface("Impinging Jet", 131) = {surf_jet[]};
Physical Surface("Radius center", 132) = {surf_axis[]};
Physical Surface("Pressure Outlet", 133) = {surf_outlet_hori[],surf_outlet_vert[]};
Physical Surface("ExtrusionPlanes", 134) = {planos_z0[], planos_zdz[]};
Physical Volume("Fluid", 135) = {todos_volumes[]};
