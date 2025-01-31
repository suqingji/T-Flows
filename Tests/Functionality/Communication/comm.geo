Rmin = 2.0;
Rmax = 2.5;
Rfilm = 2.4;
L = 16.25;

NR = 33;
N_int = NR - 4;
N_bet = 12;
N_layer_w = 9;

N_z =101;
Point(1) = {0.0, 0.0, 0.0};
Point(2) = {Rmin, 0.0, 0.0};
Point(3) = {Rfilm, 0.0, 0.0};
Point(4) = {Rmax, 0.0, 0.0};
Point(5) = {0.0, Rmin, 0.0};
Point(6) = {0.0, Rfilm, 0.0};
Point(7) = {0.0, Rmax, 0.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Circle(7) = {2, 1, 5};
Circle(8) = {3, 1, 6};
Circle(9) = {4, 1, 7};

// Interior circle
Transfinite Curve {1, 4} = N_int Using Progression 1.0;

// Intermediate circle
Transfinite Curve {2, 5} = N_bet Using Progression 0.9;

// Close to the wall
Transfinite Curve {3, 6} = N_layer_w Using Progression 0.90;

// Angular direction
Transfinite Curve {7, 8, 9} = NR;
//+
Curve Loop(1) = {1, 7, -4};
//+
Plane Surface(1) = {1};
Recombine Surface{1};

// Surface in between
Curve Loop(2) = {2, 8, -5, -7};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};
//Close to the wall
Curve Loop(3) = {3, 9, -6, -8};
Plane Surface(3) = {3};
Transfinite Surface{3};
Recombine Surface{3};

// Create Volumes
// Interior
Extrude {0, 0, L} {
  Surface {1};
  Layers{N_z};
  Recombine;
}

// Inbetween Zone
Extrude {0, 0, L} {
  Surface {2};
  Layers{N_z};
  Recombine;
}

// Close to the wall
Extrude {0, 0, L} {
  Surface {3};
  Layers{N_z};
  Recombine;
}

Physical Volume("Fluid") = {3, 2, 1};

//---------------------------
//
//  Set boundary conditions
//
//---------------------------
Physical Surface("wall") = {61};
Physical Surface("periodic-z") = {3, 2, 1, 70, 48, 26};
Physical Surface("symm-hor") = {57, 35, 17};
Physical Surface("symm-vert") = {25, 43, 65};
