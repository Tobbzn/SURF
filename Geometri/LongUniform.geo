cl1 = 0.5;
Point(1) = {0, 0, 1, cl1};
Point(2) = {0, 0, -1, cl1};
Point(3) = {10, 0, 1, cl1};
Point(4) = {-10, 0, 1, cl1};
Point(5) = {10, 0, -1, cl1};
Point(6) = {-10, 0, -1, cl1};
Circle(1) = {6, 2, 5};
Circle(2) = {5, 2, 6};
Circle(3) = {4, 1, 3};
Circle(4) = {3, 1, 4};
Line(5) = {6, 4};
Line(6) = {5, 3};
Line Loop(8) = {5, -4, -6, 2};
Ruled Surface(8) = {8};
Line Loop(10) = {5, 3, -6, -1};
Ruled Surface(10) = {10};
Line Loop(12) = {2, 1};
Plane Surface(12) = {12};
Line Loop(14) = {4, 3};
Plane Surface(14) = {14};
Surface Loop(16) = {10, 8, 14, 12};
Volume(16) = {16};