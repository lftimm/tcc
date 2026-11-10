s = 1.5;
d = 400;
l = 50;
L = d+(s*20+1)*l;
h = s*500+l;
off = 60;
half_division = off/20;

Point(1) = {0, 0, 0, 0};
Point(2) = {0, h, 0, 0};
Point(3) = {L, h, 0, 0};
Point(4) = {L, 0, 0, 0};

Point(5) = {d, 0, 0, 0};
Point(6) = {d, l, 0, 0};
Point(7) = {d+l, l, 0, 0};
Point(8) = {d+l, 0, 0, 0};

Point(9) = {d-off,0,0,0};
Point(10) = {d-off,l+off,0,0};
Point(11) = {d+l+off,l+off,0,0};
Point(12) = {d+l+off,0,0,0};

Point(13) = {d-off,h,0,0};
Point(14) = {d+l+off,h,0,0};
Point(15) = {0,l+off,0,0};
Point(16) = {L,l+off,0,0};

Point(17) = {d-half_division,0,0,0};
Point(18) = {d-half_division,l+half_division,0,0};
Point(19) = {d+l+half_division,l+half_division,0,0};
Point(20) = {d+l+half_division,0,0,0};

Point(21) = {-l/2,0,0,0};
Point(22) = {-l/2,l+off,0,0};

Line(4) = {1, 9};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(10) = {12, 11};
Line(12) = {10, 11};
Line(13) = {10, 9};
Line(15) = {4, 12};
Line(16) = {1, 15};
Line(17) = {15, 2};
Line(18) = {2, 13};
Line(19) = {13, 14};
Line(20) = {14, 3};
Line(21) = {3, 16};
Line(22) = {16, 4};
Line(23) = {16, 11};
Line(24) = {14, 11};
Line(25) = {13, 10};
Line(26) = {10, 15};
Line(27) = {5, 17};
Line(28) = {17, 9};
Line(29) = {17, 18};
Line(30) = {18, 19};
Line(31) = {19, 20};
Line(32) = {20, 8};
Line(33) = {20, 12};
Line(34) = {11, 19};
Line(35) = {18, 10};
Line(36) = {18, 6};
Line(37) = {7, 19};

Line(38) = {15, 22};
Line(39) = {1, 21};
Line(40) = {22, 21};

Curve Loop(1) = {13, -28, 29, 35};
Plane Surface(1) = {1};
Curve Loop(2) = {6, -36, -29, -27};
Plane Surface(2) = {2};
Curve Loop(3) = {7, 37, -30, 36};
Plane Surface(3) = {3};
Curve Loop(4) = {8, -32, -31, -37};
Plane Surface(4) = {4};
Curve Loop(5) = {34, 31, 33, 10};
Plane Surface(5) = {5};
Curve Loop(6) = {30, -34, -12, -35};
Plane Surface(6) = {6};
Curve Loop(7) = {23, -10, -15, -22};
Plane Surface(7) = {7};
Curve Loop(8) = {24, -23, -21, -20};
Plane Surface(8) = {8};
Curve Loop(9) = {12, -24, -19, 25};
Plane Surface(9) = {9};
Curve Loop(10) = {18, 25, 26, 17};
Plane Surface(10) = {10};
Curve Loop(11) = {4, -13, 26, -16};
Plane Surface(11) = {11};
Curve Loop(12) = {16,38,40,-39};
Plane Surface(12) = {12};

Transfinite Surface {1} = {9, 10, 18, 17};
Transfinite Surface {2} = {17, 18, 6, 5};
Transfinite Surface {3} = {6, 18, 19, 7};
Transfinite Surface {4} = {8, 7, 19, 20};
Transfinite Surface {5} = {20, 19, 11, 12};
Transfinite Surface {6} = {18, 10, 11, 19};
Transfinite Surface {7} = {12, 11, 16, 4};
Transfinite Surface {8} = {11, 14, 3, 16};
Transfinite Surface {9} = {10, 13, 14, 11};
Transfinite Surface {10} = {15, 2, 13, 10};
Transfinite Surface {11} = {1, 15, 10, 9};
Transfinite Surface {12} = {15, 22, 21, 1};


outside_house = 200;
house_eles = 100;
world_eles = 60;

a = 1.01;
b = 1;
c = 1.015;
d = 1.02;
e = 1;
f = 1;

se = 1.75;
ea = se*75;
eb = se*10;

ec = se*100;
ed = se*90;

ee = se*125;
ef = se*25;
eg = se*125;
eh = se*300;

Transfinite Curve {16, 29, 6, 8, 31, 10} = ea + 1 Using Progression f;
Transfinite Curve {13,40,22} = ea + 1 Using Progression 1/f;
Transfinite Curve {27, 36, 37, 32} = eb + 1 Using Progression 1;
Transfinite Curve {7, 30, 12, 19} = ec + 1 Using Progression 1;
Transfinite Curve {28, 35, 33} = ed + 1 Using Progression a;
Transfinite Curve {34} = ed + 1 Using Progression 1/a;
Transfinite Curve {26} = eh + 1 Using Progression b;
Transfinite Curve {4,18} = eh + 1 Using Progression 1/b;
Transfinite Curve {17} = ee + 1 Using Progression c;
Transfinite Curve {25, 24, 21} = ee + 1 Using Progression 1/c;
Transfinite Curve {20} = eg + 1 Using Progression d;
Transfinite Curve {15, 23} = eg + 1 Using Progression 1/d;
Transfinite Curve {38, 39} = ef + 1 Using Progression e;

Recombine Surface {10, 9, 11, 1, 3, 2, 6, 5, 7, 4, 8,12};

Extrude {0, 0, -1} {
  Surface{9}; Surface{10}; Surface{11};
  Surface{5}; Surface{4};  Surface{12};
  Surface{6}; Surface{1};
  Surface{2}; Surface{3};
  Surface{8}; Surface{7};
  Layers {1}; Recombine;
}

Physical Surface("Inlet", 1) = {163};
Physical Surface("Wall", 2) = {171, 93, 207, 237, 141, 123, 299};
Physical Surface("Radius", 3) = {167, 83};
Physical Surface("Slip", 4) = {71, 57, 281};
Physical Surface("Outlet", 5) = {277, 303};
Physical Surface("House", 6) = {225, 247, 137};
Physical Surface("Extrusion", 7) = {172, 12, 106, 11, 84, 10, 62, 9, 282, 8, 7, 304, 5, 128, 6, 194, 260, 3, 150, 4, 238, 2, 216, 1};
Physical Volume("Domain", 8) = {1,2,3,4,5,6,7,8,9,10,11,12};
