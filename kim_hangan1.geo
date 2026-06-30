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
Point(13) = {0.5*Djet,0,0,1.0};


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 13};
//+
Line(4) = {13, 1};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 13};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -8, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {2} = {13, 4, 5, 6};
//+
Transfinite Surface {1} = {1, 2, 3, 13};
//+

//+
Transfinite Curve {3, 8, 6} = 20 Using Bump_HWall 5;
//+
Recombine Surface {2, 1};
