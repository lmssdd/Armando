// Gmsh project created on Tue Nov 13 13:20:22 2012
x = 100.0;
y = -30.0;
z = -110.0;
h = 4.0;
Point(1) = {x + 0, y + 0, z + 0, h};
Point(2) = {x + 20, y + 0, z + 0, h};
Point(3) = {x + 40, y + 0, z + 0, h};
Point(4) = {x + 60, y + 0, z + 0, h};
Point(5) = {x + 0, y + 85, z + 0, h};
Point(6) = {x + 20, y + 85, z + 0, h};
Point(7) = {x + 40, y + 110, z + 0, h};
Point(8) = {x + 60, y + 110, z + 0, h};
BSpline(1) = {1, 2, 3, 4};
BSpline(2) = {5, 6, 7, 8};
Line(3) = {4, 8};
Line(4) = {1, 5};
Translate {-235, 0, 0} {
  Duplicata { Line{4}; }
}
Translate {-75, 0, 0} {
  Duplicata { Line{5}; }
}
/*
Translate {300, 0, 0} {
  Duplicata { Line{3}; }
}
Translate {220, 0, 0} {
  Duplicata { Line{7}; }
}
Line(9) = {11, 9};
Line(10) = {9, 1};
Line(11) = {4, 13};
Line(12) = {13, 15};
Line(13) = {16, 14};
Line(14) = {14, 8};
Line(15) = {5, 10};
Line(16) = {10, 12};
Translate {0, 50, 0} {
  Duplicata { Line{13}; }
}
Line(18) = {14, 18};
Line(19) = {16, 17};
Rotate {{0, 0, 1}, {x + 30, y + 40, z + 0}, Pi/12} {
  Line{16, 15, 4, 10, 5, 9, 6};
  Point{6, 2};
}

Translate {0, 0, 220} {
  Duplicata { Line{9, 16, 6, 5, 10, 15, 4, 1, 2, 3, 11, 14, 7, 12, 13, 18, 17, 19, 8}; }
}


Line(39) = {5, 39};
Line(40) = {56, 8};
Line Loop(41) = {2, -40, -28, -39};
Ruled Surface(42) = {41};
Line(43) = {12, 24};
Line(44) = {10, 23};
Line(45) = {11, 19};
Line(46) = {9, 20};
Line(47) = {1, 36};
Line(48) = {4, 50};
Line(49) = {13, 64};
Line(50) = {14, 67};
Line(51) = {15, 76};
Line(52) = {16, 79};
Line(53) = {18, 84};
Line(54) = {17, 87};
Line Loop(55) = {16, -6, 9, 5};
Plane Surface(56) = {55};
Line Loop(57) = {21, -22, 20, 23};
Plane Surface(58) = {57};
Line Loop(59) = {9, 46, -20, -45};
Plane Surface(60) = {59};
Line Loop(61) = {16, 43, -21, -44};
Plane Surface(62) = {61};
Line Loop(63) = {10, 4, 15, -5};
Plane Surface(64) = {63};
Line Loop(65) = {10, 47, -24, -46};
Plane Surface(66) = {65};
Line Loop(67) = {25, -23, 24, 26};
Plane Surface(68) = {67};
Line Loop(69) = {4, 2, -3, -1};
Plane Surface(70) = {69};
Line Loop(71) = {26, 28, -29, -27};
Plane Surface(72) = {71};
Line Loop(73) = {1, 48, -27, -47};
Ruled Surface(74) = {73};
Line Loop(75) = {11, 7, 14, -3};
Plane Surface(76) = {75};
Line Loop(77) = {49, -30, -48, 11};
Plane Surface(78) = {77};
Line Loop(79) = {31, -29, 30, 32};
Plane Surface(80) = {79};
Line Loop(81) = {14, -40, -31, -50};
Plane Surface(82) = {81};
Line Loop(83) = {13, -7, 12, 8};
Plane Surface(84) = {83};
Line Loop(85) = {12, 51, -33, -49};
Plane Surface(86) = {85};
Line Loop(87) = {32, -34, -38, -33};
Plane Surface(88) = {87};
Line Loop(89) = {13, 18, -17, -19};
Plane Surface(90) = {89};
Line Loop(91) = {50, 35, -53, -18};
Plane Surface(92) = {91};
Line Loop(93) = {34, 35, -36, -37};
Plane Surface(94) = {93};
Line Loop(95) = {52, 37, -54, -19};
Plane Surface(96) = {95};
Line Loop(97) = {51, 38, -52, -8};
Plane Surface(98) = {97};
*/
