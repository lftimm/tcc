dz = 0.1;
Djet = 1;
Zjet = 4*Djet;
y = 0.3*Zjet;
y2 = y;

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
Point(12) = {0.5*Djet,y,0,1.0};
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

Transfinite Curve {1, 13, 8} = 11 Using Progression 1;
Transfinite Curve {2, 7, 14} = 6 Using Progression 1;
Transfinite Curve {6, 4} = 21 Using Progression 1;
Transfinite Curve {18, 15, 3} = 9 Using Progression 1;
Transfinite Curve {17, 16, 11, 5} = 21 Using Progression 1;

Recombine Surface {5};
Recombine Surface {4};
Recombine Surface {3};
Recombine Surface {2};
Recombine Surface {1};


Extrude {0, 0, dz} {
  Surface{5}; Curve{5}; Curve{6}; Curve{11}; Surface{4}; Curve{7}; Surface{3}; Curve{16}; Curve{17}; Curve{8}; Curve{13}; Surface{2}; Curve{18}; Curve{1}; Surface{1}; Curve{2}; Curve{14}; Curve{3}; Curve{15}; Curve{4}; Layers {1}; Recombine;
}
