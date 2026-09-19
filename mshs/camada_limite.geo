d = 400;
l = 50;
L = d+21*l;
h = 500+l;
off = 60;

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

Line(4) = {1, 9};
Line(5) = {9, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 12};
Line(10) = {12, 11};
Line(11) = {11,7};
Line(12) = {10, 11};
Line(13) = {10, 9};
Line(14) = {6, 10};
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




Curve Loop(1) = {12, 11, -7, 14};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 14, 13};
Plane Surface(2) = {2};

Curve Loop(3) = {11, 8, 9, 10};
Plane Surface(3) = {3};

Curve Loop(4) = {26, -16, 4, -13};
Plane Surface(4) = {4};
Curve Loop(5) = {25, 26, 17, 18};
Plane Surface(5) = {5};
Curve Loop(6) = {19, 24, -12, -25};
Plane Surface(6) = {6};
Curve Loop(7) = {20, 21, 23, -24};
Plane Surface(7) = {7};
Curve Loop(8) = {15, 10, -23, 22};
Plane Surface(8) = {8};

Transfinite Surface {1} = {10, 11, 7, 6};
Transfinite Surface {2} = {9, 10, 6, 5};
Transfinite Surface {3} = {7, 11, 12, 8};
Transfinite Surface {4} = {1, 15, 10, 9};
Transfinite Surface {5} = {15, 2, 13, 10};
Transfinite Surface {6} = {10, 13, 14, 11};
Transfinite Surface {7} = {11, 14, 3, 16};
Transfinite Surface {8} = {12, 11, 16, 4};


house_eles = 200;
prog = 1.003;
Transfinite Curve {14, 9} = house_eles Using Progression prog;
Transfinite Curve {5,11} = house_eles Using Progression 1/prog;

outside_house = 300;
Transfinite Curve {13, 16, 6, 8, 10, 22} = outside_house Using Progression 1;

Transfinite Curve {7, 12, 19} = outside_house Using Progression 1;

world_eles = 100;
world_prog = 1.02;
Transfinite Curve {4,  18} = world_eles Using Progression world_prog;

Transfinite Curve {15, 23} = world_eles Using Progression world_prog;

Transfinite Curve {20} = world_eles Using Progression 1/world_prog;

Transfinite Curve {21, 24, 25} = world_eles Using Progression 1/world_prog;

Transfinite Curve {17, 26} = world_eles Using Progression world_prog;


Recombine Surface {4, 2, 1, 8, 3, 7, 6, 5};