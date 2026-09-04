// -----------------------------
// AIRFOIL CFD-ready 2D mesh
// Re ~ 6.7e05 (U = 10 m/s)
// -----------------------------

// SetFactory("OpenCASCADE");
SetFactory("Built-in");

// ----- Global Settings -----
General.NumThreads = 0;

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;

Mesh.Smoothing = 10;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
// ---------------------------

// NOTE: these values don't matter if defined already in meshEEW.sh. These are merely meant as a failsafe
// ----- Parameters -----
True = 1;
False = 0;

If (!Exists(AoA))
    AoA = 10;
EndIf

If (!Exists(chord))
    chord = 1.0;
EndIf

If (!Exists(span))
    span = 5;
EndIf

If (!Exists(wingspan))
    wingspan = 5;
EndIf

If (!Exists(far))
    far = 5;
EndIf

If (!Exists(firstlayer))
    firstlayer = 1e-5;
EndIf

If (!Exists(layers))
    layers = 300; //default layers=300
EndIf

If (!Exists(wing_layers))
    wing_layers = 200; //default layers=200
EndIf

If (!Exists(front))
    front = far;
EndIf

If (!Exists(rear))
    rear = 3*far;
EndIf

If (!Exists(dim))
    dim = 2; // 2 or 3
EndIf

If (!Exists(filled))
    filled = False; // 0(False) or 1(True)
EndIf

If (!Exists(quads))
    quads = False; // 0(False) or 1(True)
EndIf

If (!Exists(groundres))
    groundres = 200;
EndIf

If (!Exists(wakeres))
    wakeres = 300;
EndIf

If (!Exists(inletres))
    inletres = 40;
EndIf

If (!Exists(outletres))
    outletres = 20;
EndIf

// ----------------------


AoA_rad = AoA * Pi/180;
Printf("AoA(deg) = %g AoA(rad) = %g firstlayer = %g", AoA,AoA_rad,firstlayer);

// Global mesh controls
Mesh.CharacteristicLengthMin = 1e-6;
Mesh.CharacteristicLengthMax = 0.5;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.MinimumElementsPerTwoPi = 50;


//----- Import Airfoil -----
Merge "Airfoil_points.geo";

If (!Exists(BLAirfoil))
    BLAirfoil = 0.1 * chord;
EndIf
If (!Exists(BLGround))
    BLGround = 0;
EndIf
If (!Exists(hc)) // this represents h/c for the ground height
    hc = 1;
EndIf
ground = hc * chord;


// ----- SCALE (chord) -----
Dilate {{0,0,0}, {chord, chord, 1}} {
    Curve{1,2,3,4};
}

// ----- ROTATE (AoA) -----
//Rotate {{0,0,1}, {0,0,0}, -AoA_rad} {
//    Curve{1,2,3,4};
//}

OutUpper = 30000;
WakeUpper = TEoff_up;
WakeLower = TEoff_low;
OutLower = 30006;
OutMid = 30007;


tUp = BLAirfoilUp * chord;
tLow = BLAirfoilLow * chord;
// Origin shifted to TE(look at dat_to_geo for more info)!!
xTE = 0;
yTE = 0;

Point(OutUpper) = {rear,yTE + tUp,0};
Point(OutLower) = {rear,yTE - tLow,0};
Point(OutMid) = {rear,yTE,0};


//Lines starting from 50 in the structured part
Line(50) = {OutUpper,OutMid};
Line(51) = {WakeUpper,OutUpper};
Line(56) = {WakeLower,OutLower};
Line(60) = {LEoff,LEpoint};
Line(61) = {OutMid,TEpoint};
Line(62) = {TEpoint,WakeUpper};
Line(63) = {TEpoint,WakeLower};
Line(57) = {OutLower,OutMid};

//--------------------------



// ----- FARFIELD -----
Point(20001) = {rear, far,0};
Point(20002) = {-front, far,0};
Point(20003) = {-front,yTE-ground,0};
Point(20004) = {rear,yTE-ground,0};

//Point(20003) = {-front,yTE-ground+BLGround,0};
//Point(20004) = {rear,yTE-ground+BLGround,0};
//Point(20013) = {-front,yTE-ground,0};
//Point(20014) = {rear,yTE-ground,0};

//this is inside Airfoil_points.geo 
//assumes standard airfoil dat file structue TE-LE-TE
Physical Point("TE") = {TEpoint};
Physical Point("LE") = {LEpoint};

Line(21) = {20001,20002}; //farfield
Line(22) = {20002,20003}; //inlet
//Lines from the structured part
Line(25) = {OutUpper,20001};

Transfinite Curve{22} = inletres;
Transfinite Curve{25} = outletres; // NOTE: this is only for part of the outlet


//ground
//Line(213) = {20004,20014};
//Line(223) = {20014,20013};
//Line(233) = {20013,20003};

If (BLAirfoilUp>=ground)
    Line(23) = {20003,GNDpoint};
    Transfinite Curve{23} = Floor(groundres/2) Using Progression 0.98;
    Curve Loop(100) = {21,22,23,-4,-3,51,25}; //Unstructured outer domain block
    outlet_curves[]   = {57,-50,25};
    ground_curves[]   = {23,5,56};
Else
    Line(23) = {20003,20004}; //ground-BLGround
    Line(24) = {20004,OutLower};
    Transfinite Curve{23} = groundres;
    Curve Loop(100) = {21,22,23,24,-56,-5,-4,-3,51,25}; //Unstructured outer domain block
    outlet_curves[]   = {24,57,-50,25};
    ground_curves[]   = {23};
EndIf

Compound Curve{4,5};
// airfoil edges are 1 & 2
Curve Loop(200) = {60,-1,62,3}; //BLupper domain block
Curve Loop(300) = {60,2,63,-5,-4}; //BLlower domain block
Curve Loop(400) = {61,62,51,50}; //upperwake domain block
Curve Loop(500) = {61,63,56,57}; //lowerwake domain block

Curve Loop(999) = {1,2}; //airfoil block


// Fluid surface with hole
Plane Surface(30) = {100}; //outer surface
Plane Surface(40) = {200}; //BLupper surface
Plane Surface(41) = {300}; //BLlower surface
Plane Surface(42) = {400}; //upperwake surface
Plane Surface(43) = {500}; //lowerwake surface
Plane Surface(99) = {999}; //airfoil surface

AirfoilBump = 0.25;
Transfinite Curve{-1} = N_UP Using Bump AirfoilBump;
Transfinite Curve{-3} = N_UP Using Progression 0.99;
Transfinite Curve{2} = N_LOW Using Bump AirfoilBump;
Transfinite Curve{4} = N4; // Using Progression 0.97;
Transfinite Curve{5} = N5 Using Progression 0.98;
Transfinite Curve{51,56,-61} = wakeres Using Progression 1.01;
Transfinite Curve{51} = wakeres Using Progression 1.005;
Transfinite Curve{63,-57} = NUMlayers Using Progression grTElow;
Transfinite Curve{62,-50} = NUMlayers Using Progression grTEup;
Transfinite Curve{60} = NUMlayers Using Progression 2-grLE;

//Transfinite Curve{23,223} = 200;
//Transfinite Curve{213,233} = 50;

Transfinite Surface(40); //BL surface
Transfinite Surface(41) = {LEoff, LEpoint, TEpoint, WakeLower}; //BL surface
Transfinite Surface(42); //upperwake surface
Transfinite Surface(43); //lowerwake surface

// Recombining to make Quads, if needed
If (quads == True)
    Recombine Surface {40,41,42,43};
EndIf


//______ LOWER FIELDS _________


// Farfield base size
Field[3] = MathEval;
Field[3].F = "0.2";

// Combine fields
Field[4] = Min;
Field[4].FieldsList = {3};

Background Field = 4;

//_____________________________

If ((dim == 3) && (filled == True))

    // ----- EXTRUDE -----

    If (wingspan == span)
        
        Extrude {0,0,span*chord}
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            Surface{99};
            Layers{layers};
            Recombine;
        }

    Else 
    
        Extrude {0,0,-(span-wingspan)*chord} //NOTE: span > wingspan always
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            Surface{99};
            Layers{layers-wing_layers};
            Recombine;
        }
   
        Extrude {0,0,wingspan*chord}
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            //Surface{99};
            Layers{wing_layers};
            Recombine;
        }
        
    EndIf

    // ----- PHYSICAL GROUPS -----

    // ----- Find external surfaces -----

    tol = 1e-8 * chord;

    // x = -front
    inlet[] =
        Surface In BoundingBox{
            -front-tol, -far-tol, -(span-wingspan)-tol,
            -front+tol,  far+tol, wingspan+tol
        };

    // x = rear
    outlet[] =
        Surface In BoundingBox{
            rear-tol, yTE-ground-tol, -(span-wingspan)-tol,
            rear+tol, far+tol, wingspan+tol
        };

    // y = far
    topwall[] =
        Surface In BoundingBox{
            -front-tol, far-tol, -(span-wingspan)-tol,
            rear+tol,  far+tol, wingspan+tol
        };

    // y = yTE-ground
    ground[] =
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, -(span-wingspan)-tol,
            rear+tol,  yTE-ground+tol, wingspan+tol
        };

    // z = 0
    backSurf[] =
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, -(span-wingspan)-tol,
            rear+tol,  far+tol, -(span-wingspan)+tol
        };

    // z = span
    frontSurf[] = //because extrusion is done along positive axis!
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, wingspan-tol,
            rear+tol,  far+tol, wingspan+tol
        };
        
    // ----- PHYSICAL GROUPS -----

    Physical Volume("fluid") = {1,2,3,4,5};

    Physical Surface("inlet")   = inlet[];
    Physical Surface("outlet")  = outlet[];
    Physical Surface("topwall") = topwall[];
    Physical Surface("ground")  = ground[];
    
    // These IDs correspond to the airfoil side surfaces created by the extrusion.
    // If the extrusion topology changes, update them using Gmsh's Visibility tool.
    // Physical Surface("airfoil")   = {1059, 1081};

    Physical Surface("front")   = frontSurf[];
    Physical Surface("back")    = backSurf[];



ElseIf ((dim == 3) && (filled == False))

    // ----- EXTRUDE -----

    If (wingspan == span)
        
        Extrude {0,0,-span*chord}
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            //Surface{99};
            Layers{layers};
            Recombine;
        }

    Else 
    
        Extrude {0,0,(span-wingspan)*chord}
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            Surface{99};
            Layers{layers-wing_layers};
            Recombine;
        }
   
        Extrude {0,0,-wingspan*chord}
        {
            Surface{30};
            Surface{40};
            Surface{41};
            Surface{42};
            Surface{43};
            //Surface{99};
            Layers{wing_layers};
            Recombine;
        }
        
    EndIf

    // ----- PHYSICAL GROUPS -----

    // ----- Find external surfaces -----

    tol = 1e-8 * chord;

    // x = -front
    inlet[] =
        Surface In BoundingBox{
            -front-tol, -far-tol, -wingspan-tol,
            -front+tol,  far+tol, (span-wingspan)+tol
        };

    // x = rear
    outlet[] =
        Surface In BoundingBox{
            rear-tol, yTE-ground-tol, -wingspan-tol,
            rear+tol, far+tol, (span-wingspan)+tol
        };

    // y = far
    topwall[] =
        Surface In BoundingBox{
            -front-tol, far-tol, -wingspan-tol,
            rear+tol,  far+tol, (span-wingspan)+tol
        };

    // y = yTE-ground
    ground[] =
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, -wingspan-tol,
            rear+tol,  yTE-ground+tol, (span-wingspan)+tol
        };

    // z = 0
    backSurf[] =
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, -wingspan-tol,
            rear+tol,  far+tol, -wingspan+tol
        };

    // z = span
    frontSurf[] = //because extrusion is done along positive axis!
        Surface In BoundingBox{
            -front-tol, yTE-ground-tol, (span-wingspan)-tol,
            rear+tol,  far+tol, (span-wingspan)+tol
        };
        
    // ----- PHYSICAL GROUPS -----

    Physical Volume("fluid") = {1,2,3,4,5};

    Physical Surface("inlet")   = inlet[];
    Physical Surface("outlet")  = outlet[];
    Physical Surface("topwall") = topwall[];
    Physical Surface("ground")  = ground[];

    // These IDs correspond to the airfoil side surfaces created by the extrusion.
    // If the extrusion topology changes, update them using Gmsh's Visibility tool.
    Physical Surface("airfoil")   = {1059, 1081};

    Physical Surface("front")   = frontSurf[];
    Physical Surface("back")    = backSurf[];


// ----------------------

ElseIf ((dim == 2) && (filled == True))

    // ----- PHYSICAL GROUPS -----

    Color {0,220,220} { Surface{30}; }   // outer
    Color {255,150,150} { Surface{40}; }   // BL upper
    Color {150,150,255} { Surface{41}; }   // BL lower
    Color {150,255,150} { Surface{42}; }   // wake upper
    Color {255,255,150} { Surface{43}; }   // wake lower

    Physical Surface("fluid") = {30,40,41,42,43,99};
    Physical Curve("outlet")   = outlet_curves[];
    Physical Curve("ground")   = ground_curves[];
    Physical Curve("inlet")    = {22};
    Physical Curve("topwall")  = {21};
    Physical Curve("airfoil") = {1,2};

ElseIf ((dim == 2) && (filled == False))

    // ----- PHYSICAL GROUPS -----

    Color {0,220,220} { Surface{30}; }   // outer
    Color {255,150,150} { Surface{40}; }   // BL upper
    Color {150,150,255} { Surface{41}; }   // BL lower
    Color {150,255,150} { Surface{42}; }   // wake upper
    Color {255,255,150} { Surface{43}; }   // wake lower

    Physical Surface("fluid") = {30,40,41,42,43};
    Physical Curve("outlet")   = outlet_curves[];
    Physical Curve("ground")   = ground_curves[];
    Physical Curve("inlet")    = {22};
    Physical Curve("topwall")  = {21};
    Physical Curve("airfoil") = {1,2};

Else

    Error("Invalid dim/filled value");

EndIf
