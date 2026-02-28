
Lc = 0.01;

p1 = newp; Point(p1) = {0.0, 0.0, 0.0, Lc};
p2 = newp; Point(p2) = {1.0, 0.0, 0.0, Lc};
p3 = newp; Point(p3) = {1.0, 0.0, 1.0, Lc};
p4 = newp; Point(p4) = {0.0, 0.0, 1.0, Lc};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

ll = newll; Line Loop(ll) = {l1,l2,l3,l4};

ps = news;  Plane Surface(ps) = {ll};

Ext[] = Extrude{0, Lc, 0}{Surface{ps}; Layers{1}; Recombine;};

Physical Surface("left")   = {27};
Physical Surface("right")  = {19};
Physical Surface("top")    = {23};
Physical Surface("bottom") = {15};
Physical Surface("front")  = {6};
Physical Surface("back")   = {28};

Physical Volume("all")     = {1};

Mesh.Format = 10;
Mesh.CharacteristicLengthFactor = 1;
Mesh.CharacteristicLengthMin = 0.0;
Mesh.CharacteristicLengthMax = 1.0e22;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 1;

Mesh.ElementOrder = 1;
Mesh.Explode = 1.0;
//Mesh.SecondOrderIncomplete = 1;
Mesh.SaveGroupsOfNodes = 1;

Mesh.Algorithm = 5;
// 1: MeshAdapt, 2:Automatic, 5:Delaunay, 6:Frontal, 7:bamg, 8:DelQuad, 9: Packing of parallelograms

Mesh.RandomFactor = 1e-9;
Mesh.RecombinationAlgorithm = 1;		// 0: standard, 1:blossom
Mesh.RecombineAll = 1;
Mesh.Smoothing = 3;
Mesh.Optimize = 1;
Mesh.RefineSteps = 5;
