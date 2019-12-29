//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0.5, 0, 1.0};
//+
Point(3) = {0.5, 0.5, 0, 1.0};
//+
Point(4) = {0.5, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wall") = {4, 2};
//+
Physical Curve("in", 6) = {1};
//+
Physical Curve("out", 7) = {3};
//+
Physical Surface("plane", 8) = {1};
