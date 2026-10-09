d = 400;
l = 50;
L = d+21*l;
h = 500+l;
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

Physical Surface("Wall", 1) = {200, 134, 222, 182, 192, 274, 120};
Physical Surface("Inlet", 2) = {102,80};
Physical Surface("Outlet", 3) = {278, 252};
Physical Surface("Slip", 4) = {68,54,256};
Physical Surface("Extrude", 5) = {103,81,59,191,169,125,147,235,213,257,279,1,2,3,4,5,6,7,8,9,10,11};
Physical Volume("Fluid",6) = {1,2,3,4,5,6,7,8,9,10,11};

outside_house = 200;
house_eles = 100;
world_eles = 60;

a = 1.015;
b = 1.02;
c = 1.02;
d = 1.03;

ea = 75;
eb = 6;

ec = 100;
ed = 50;

ee = 100;

Transfinite Curve {16, 13, 29, 6, 8, 31, 10, 22} = ea Using Progression 1;
Transfinite Curve {27, 36, 37, 32} = eb Using Progression 1;
Transfinite Curve {7, 30, 12, 19} = ec Using Progression 1;

Transfinite Curve {28, 35, 33} = ed Using Progression a;
Transfinite Curve {34} = ed Using Progression 1/a;

Transfinite Curve {26} = ee Using Progression b;
Transfinite Curve {4,18} = ee Using Progression 1/b;

Transfinite Curve {17} = ee Using Progression c;
Transfinite Curve {25, 24, 21} = ee Using Progression 1/c;

Transfinite Curve {20} = ee Using Progression d;
Transfinite Curve {15, 23} = ee Using Progression 1/d;

Recombine Surface {10, 9, 11, 1, 3, 2, 6, 5, 7, 4, 8};

Extrude {0, 0, 1} {
  Surface{9}; Surface{10}; Surface{11};
  Surface{5}; Surface{4};
  Surface{6}; Surface{1};
  Surface{2}; Surface{3};
  Surface{8}; Surface{7};
  Layers {1}; Recombine;
}

