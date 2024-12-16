//+
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = 0.25;
// Control points
// NOTE: See the cp_vane at the bottom of the file
cp_wall = 0.5;
cp_vol = 0.2;

Mesh.RandomFactor = 5.e-10;
// Choose 2D Algorithm
// 1 : MeshAdapt
// 2 : Automatic
// 5 : Delaunay
// 6 : Frontal-Delaunay
Mesh.Algorithm =1;
//Choose 3D Algorithm
// 1 : Delaunay
Mesh.Algorithm3D =1;
//
vane_diameter = 0.56418958354; // Pi*r^2 = 1 => r = Sqrt[(1/Pi)]
//radius = 4 * vane_diameter; // 2 * 2*r
radius = 2.5;
xmid = 0.75 * vane_diameter * 2;
ymid = 0;
zmid = 0.75 * vane_diameter * 2;
vane_start = 0.25;//xmid - (vane_diameter/2)*2;
Printf("%g", vane_start);
Printf("%g", xmid);
Printf("%g", vane_start+vane_diameter*2);

vane_thickness = 0.05;  // total_vane_thickness 0.005*2
lx = vane_start + vane_diameter; 
ly = vane_diameter; 
lz = vane_start + vane_diameter;

//+
// RIGHT Vane
// Draw circle on hot side and convert it to surface
c1 = newc; Circle(c1) = {lx, 0, vane_thickness, vane_diameter, 0, 2*Pi};
// Create Hot surface to apply boundary condition
rh_cl = newcl; Curve Loop(rh_cl) = {c1};
rh_s = news;  Plane Surface(rh_s) = {rh_cl};
// Extrude surface from Hot to cold
right_vane_extrude[] = Extrude {0, 0, -2*vane_thickness} { Surface{rh_s}; };
// Create Cold surface to apply boundary condition
rc_s= news; rc_s = {right_vane_extrude[0]};
// Create Transition surface from Hot to Cold 
rt_s = news; rt_s = {right_vane_extrude[2]};
// Define the right vane volume
r_v = newv; r_v = {right_vane_extrude[1]};
//
//Printf("%g", r_v);
//Printf("%g", #right_vane_extrude() );
//Printf("%g %g %g", right_vane_extrude[0], right_vane_extrude[1],right_vane_extrude[2]);


//DOWN Vane
// Draw circle on hot side
c2 = newc; Circle(c2) = {vane_thickness, 0, -lz, vane_diameter, 0, 2*Pi};
Rotate {{0, 1, 0}, {vane_thickness, 0, -lz}, Pi/2} { Curve{c2}; }
// Create Hot surface to apply boundary condition
dh_cl = newcl; Curve Loop(dh_cl) = {c2}; 
dh_s = news; Plane Surface(dh_s) = {dh_cl};
// Extrude circle from Hot to cold
down_vane_extrude[] = Extrude {-2*vane_thickness, 0, 0 } {Surface{dh_s};};
// Create Cold surface to apply boundary condition
dc_s = news; dc_s = {down_vane_extrude[0]};
// Create Transition surface from Hot to Cold 
dt_s = news; dt_s = {down_vane_extrude[2]};
// Define the down vane volume
d_v = newv; d_v = {down_vane_extrude[1]};
Printf("%g %g %g", down_vane_extrude[0], down_vane_extrude[1],down_vane_extrude[2]);
//
// LEFT vane
// Draw circle on hot side
c3 = newc; Circle(c3) = {-lx, 0, -vane_thickness, vane_diameter, 0, 2*Pi};
// Create Hot surface to apply boundary condition
lh_cl = newcl; Curve Loop(lh_cl) = {c3}; 
lh_s = news; Plane Surface(lh_s) = {lh_cl};
// Extrude Surface from Hot to cold
left_vane_extrude[] = Extrude {0,0,+2*vane_thickness}{Surface{lh_s};};
// Create Cold surface to apply boundary condition
lc_s = news; lc_s = {left_vane_extrude[0]};
// Create Transition surface from Hot to Cold 
lt_s = news; lt_s = {left_vane_extrude[2]};
// Define the left vane volume
l_v = newv; l_v = {left_vane_extrude[1]};

// UP vane
// Draw circle on hot side
c4 = newc; Circle(c4) = {-vane_thickness, 0, lz, vane_diameter, 0, 2*Pi};
Rotate {{0, 1, 0}, {-vane_thickness, 0, lz}, Pi/2} { Curve{c4}; }
// Create Hot surface to apply boundary condition
uh_cl = newcl; Line Loop(uh_cl) = {c4}; 
uh_s = news; Plane Surface(uh_s) = {uh_cl};
// Extrude circle from Hot to cold
up_vane_extrude[] = Extrude {2*vane_thickness, 0, 0 } { Surface{uh_s}; };
// Create Cold surface to apply boundary condition
uc_s = news; uc_s = {up_vane_extrude[0]};
// Create Transition surface from Hot to Cold 
ut_s = news; ut_s = {up_vane_extrude[2]};
// Create Up vane volume
u_v = newv; u_v = {up_vane_extrude[1]};
//
//Physical Surfaces
//Right vane
Physical Surface("RightHot",188) = {rh_s};
Physical Surface("RightCold",183) = {rc_s};
Physical Surface("RightTransition",1820) = {rt_s};
// DOWN VANE
Physical Surface("DownHot", 48) = {dh_s};
Physical Surface("DownCold",43) = {dc_s};
Physical Surface("DownTransition",420) = {dt_s};
//// LEFT VANE
Physical Surface("LeftHot",128) = {lh_s};
Physical Surface("LeftCold",123) = {lc_s};
Physical Surface("LeftTransition",1220) = {lt_s};
// UP VANE
Physical Surface("UpHot", 218) = {uh_s};
Physical Surface("UpCold",213) = {uc_s};
Physical Surface("UpTransition",2120) = {ut_s};
//
////
//r_sl = newsl; Surface Loop(r_sl) = {rh_s, rt_s, rc_s};
//d_sl = newsl; Surface Loop(d_sl) = {dh_s, dt_s, dc_s};
//l_sl = newsl; Surface Loop(l_sl) = {lh_s, lt_s, lc_s};
//u_sl = newsl; Surface Loop(u_sl) = {uh_s, ut_s, uc_s};
////
//r_v = newv; Volume(r_v) = {r_sl};
//d_v = newv; Volume(d_v) = {d_sl};
//l_v = newv; Volume(l_v) = {l_sl};
//u_v = newv; Volume(u_v) = {u_sl};
////+
glass= newv; Sphere(glass) = {0, 0, 0, radius};
//// Introduce new points on the sphere to control mesh size- coarse mesh
cpr = newp; Point(cpr) = {radius, 0, 0, cp_wall};
cpl = newp; Point(cpl) = {-radius, 0,0, cp_wall};
cpd = newp; Point(cpd) = {0, 0,-radius, cp_wall};
cpu = newp; Point(cpu) = {0, 0, radius, cp_wall};
cpt = newp; Point(cpt) = {0, radius,0, cp_wall};
cpb = newp; Point(cpb) = {0,-radius,0, cp_wall};
Point{cpr} In Volume {glass};
Point{cpl} In Volume {glass};
Point{cpu} In Volume {glass};
Point{cpd} In Volume {glass};
Point{cpt} In Volume {glass};
Point{cpb} In Volume {glass};
//
cpn1 = newp; Point(cpn1) = {0.85*radius, 0, 0, cp_vol};
cpn2 = newp; Point(cpn2) = {-0.85*radius, 0, 0, cp_vol};
cpn3 = newp; Point(cpn3) = {0, 0.85*radius, 0, cp_vol};
cpn4 = newp; Point(cpn4) = {0, -0.85*radius, 0, cp_vol};
cpn5 = newp; Point(cpn5) = {0, 0, 0.85*radius, cp_vol};
cpn6 = newp; Point(cpn6) = {0, 0, -0.85*radius, cp_vol};
cpn7 = newp; Point(cpn7) = {xmid, 0.65*radius, 0, cp_vol};
cpn8 = newp; Point(cpn8) = {-xmid, 0.65*radius, 0, cp_vol};
cpn9 = newp; Point(cpn9) = {xmid, -0.65*radius, 0, cp_vol};
cpn10= newp; Point(cpn10) = {-xmid, -0.65*radius, 0, cp_vol};
cpn11 = newp; Point(cpn11) = {xmid, 0, 0.65*radius, cp_vol};
cpn12 = newp; Point(cpn12) = {xmid, 0, -0.65*radius, cp_vol};
cpn13 = newp; Point(cpn13) = {-xmid, 0, 0.65*radius, cp_vol};
cpn14 = newp; Point(cpn14) = {-xmid, 0, -0.65*radius, cp_vol};
Point{cpn1} In Volume {glass};
Point{cpn2} In Volume {glass};
Point{cpn3} In Volume {glass};
Point{cpn4} In Volume {glass};
Point{cpn5} In Volume {glass};
Point{cpn6} In Volume {glass};
Point{cpn7} In Volume {glass};
Point{cpn8} In Volume {glass};
Point{cpn9} In Volume {glass};
Point{cpn10} In Volume {glass};
Point{cpn11} In Volume {glass};
Point{cpn12} In Volume {glass};
Point{cpn13} In Volume {glass};
Point{cpn14} In Volume {glass};

//
Physical Surface("Wall", 23) = {21};
//////
Physical Volume ("VOLUME",8888) = { BooleanDifference{ Volume{glass}; Delete; }{ Volume{r_v}; Volume{d_v}; Volume{l_v}; Volume{u_v};} };
////
//Physical Volume ("VOLUME",8888) = {BooleanDifference{ Volume{glass};
//Delete; }{Volume{r_v};}};
Characteristic Length { 1,2,3,4,5,6,7,8,9,10} = 0.02;
