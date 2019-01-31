// Gmsh project created on Tue Nov 13 13:20:22 2012
x = 0.000;
y = 0.000;
z = -0.110;
h = 0.03;
Point(1) = {x + 0.000, y + 0.000, z + 0.000, h};
Point(2) = {x + 0.020, y + 0.000, z + 0.000, h};
Point(3) = {x + 0.040, y + 0.000, z + 0.000, h};
Point(4) = {x + 0.060, y + 0.000, z + 0.000, h};
Point(5) = {x + 0.000, y + 0.085, z + 0.000, h};
Point(6) = {x + 0.020, y + 0.085, z + 0.000, h};
Point(7) = {x + 0.040, y + 0.110, z + 0.000, h};
Point(8) = {x + 0.060, y + 0.110, z + 0.000, h};
BSpline(1) = {1, 2, 3, 4};
BSpline(2) = {5, 6, 7, 8};
Line(3) = {4, 8};
Line(4) = {1, 5};
Translate {-0.235, 0.000, 0.000} {
  Duplicata { Line{4}; }
}
Line(6) = {9, 1};
Line(7) = {10, 5};
Translate {-0.1, 0, 0.0} {
  Duplicata { Line{5}; }
}
Line(9) = {12, 10};
Line(10) = {11, 9};
Point(13) = {x-0.335, y-0.25, z, 1.0};
Rotate {{0, 0, 1}, {x-0.335, y-0.25, 0}, Pi*1/12} {
  Duplicata { Line{8}; }
}
Circle(12) = {11, 13, 14};
Circle(13) = {12, 13, 15};
Rotate {{0, 0, 1}, {x+0.03, y, z}, Pi/12} {
  Line{4, 7, 6, 5, 9, 10, 8, 13, 12, 11};
  Point{2, 6};
}
Translate {0, 0, 0.220} {
  Duplicata { Line{11, 13, 8, 12, 10, 9, 5, 6, 7, 4, 1, 2, 3}; }
}
Line(27) = {14, 16};
Line(28) = {11, 25};
Line(29) = {9, 35};
Line(30) = {1, 47};
Line(31) = {4, 61};
Line(32) = {15, 17};
Line(33) = {12, 20};
Line(34) = {10, 39};
Line(35) = {5, 51};
Line(36) = {8, 67};
Line Loop(37) = {12, 27, -17, -28};
Ruled Surface(38) = {37};
Line Loop(39) = {10, 29, -18, -28};
Ruled Surface(40) = {39};
Line Loop(41) = {6, 30, -21, -29};
Ruled Surface(42) = {41};
Line Loop(43) = {1, 31, -24, -30};
Ruled Surface(44) = {43};
Line Loop(45) = {13, 32, -15, -33};
Ruled Surface(46) = {45};
Line Loop(47) = {9, 34, -19, -33};
Ruled Surface(48) = {47};
Line Loop(49) = {2, 36, -25, -35};
Ruled Surface(50) = {49};
Line Loop(51) = {11, -13, -8, 12};
Ruled Surface(52) = {51};
Line Loop(53) = {14, -15, -16, 17};
Ruled Surface(54) = {53};
Line Loop(55) = {8, 9, -5, -10};
Ruled Surface(56) = {55};
Line Loop(57) = {16, 19, -20, -18};
Ruled Surface(58) = {57};
Line Loop(59) = {5, 7, -4, -6};
Ruled Surface(60) = {59};
Line Loop(61) = {20, 22, -23, -21};
Ruled Surface(62) = {61};
Line Loop(63) = {4, 2, -3, -1};
Ruled Surface(64) = {63};
Line Loop(65) = {23, 25, -26, -24};
Ruled Surface(66) = {65};
Extrude {1.0, 0, 0} {
  Line{31, 26, 36, 3};
}
Extrude {{0, 0, 1}, {x+1.06, y+0.4, z}, Pi/2} {
  Line{79, 67, 71, 75};
}
