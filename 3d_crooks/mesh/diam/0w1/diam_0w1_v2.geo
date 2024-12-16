//+
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = 0.1;
cp_vane = 0.01;
Mesh.RandomFactor = 5.e-10;
// Choose 2D Algorithm
// 1 : MeshAdapt
// 2 : Automatic
// 5 : Delaunay
// 6 : Frontal-Delaunay
Mesh.Algorithm = 1;
// Choose 3D Algorithm
// 1: Delaunay 
Mesh.Algorithm3D = 1; 
//
vane_diagonal = 1.414235;
//radius = 2 * vane_diagonal;
radius = 2.5;
//xmid = 0.75 * vane_diagonal;
ymid = 0;
//zmid = 0.75 * vane_diagonal;
vane_start = 0.25; //xmid - (vane_diagonal/2);
xmid = vane_start + (vane_diagonal/2);
zmid = vane_start + (vane_diagonal/2);
vane_thickness = 0.05; 
lx = vane_start + vane_diagonal; 
ly = vane_diagonal; 
lz = vane_start + vane_diagonal;

// Right Vane
// (Hot) POINTS

p01 = newp; Point(p01) = {+vane_start, +ymid, +vane_thickness, cp_vane};
p02 = newp; Point(p02) = {+xmid, +ly/2, +vane_thickness, cp_vane};
p03 = newp; Point(p03) = {+lx, +ymid, +vane_thickness, cp_vane};
p04 = newp; Point(p04) = {+xmid, -ly/2, +vane_thickness, cp_vane};
// (Cold) Points
p05 = newp; Point(p05) = {+vane_start, +ymid, -vane_thickness, cp_vane};
p06 = newp; Point(p06) = {+xmid, +ly/2, -vane_thickness, cp_vane};
p07 = newp; Point(p07) = {+lx, +ymid, -vane_thickness, cp_vane};
p08 = newp; Point(p08) = {+xmid, -ly/2, -vane_thickness, cp_vane};
// Bottom Vane
// (Hot) POINTS
p09 = newp; Point(p09) = { +vane_thickness,+ymid, -vane_start, cp_vane};
p10 = newp; Point(p10) = { +vane_thickness,+ly/2, -zmid, cp_vane};
p11 = newp; Point(p11) = { +vane_thickness,+ymid, -lz, cp_vane};
p12 = newp; Point(p12) = { +vane_thickness,-ly/2,-zmid, cp_vane};
// (Cold) Points
p13 = newp; Point(p13) = { -vane_thickness, +ymid, -vane_start, cp_vane};
p14 = newp; Point(p14) = { -vane_thickness, +ly/2, -zmid,cp_vane};
p15 = newp; Point(p15) = { -vane_thickness, +ymid, -lz, cp_vane};
p16 = newp; Point(p16) = { -vane_thickness, -ly/2, -zmid, cp_vane};
// Left Vane
//(HOT) Points
p17 = newp; Point(p17) = {-vane_start, +ymid, -vane_thickness, cp_vane};
p18 = newp; Point(p18) = {-xmid, +ly/2, -vane_thickness, cp_vane};
p19 = newp; Point(p19) = {-lx, +ymid, -vane_thickness, cp_vane};
p20 = newp; Point(p20) = {-xmid, -ly/2, -vane_thickness, cp_vane};
//(COLD) POINTS
p21 = newp; Point(p21) = {-vane_start, +ymid, +vane_thickness, cp_vane};
p22 = newp; Point(p22) = {-xmid, +ly/2, +vane_thickness, cp_vane};
p23 = newp; Point(p23) = {-lx, +ymid, +vane_thickness, cp_vane};
p24 = newp; Point(p24) = {-xmid, -ly/2, +vane_thickness, cp_vane};
// Top Vane
// (HOT) Points
p25 = newp; Point(p25) = { -vane_thickness, +ymid, vane_start, cp_vane};
p26 = newp; Point(p26) = { -vane_thickness, +ly/2, zmid,cp_vane};
p27 = newp; Point(p27) = { -vane_thickness, +ymid, lz, cp_vane};
p28 = newp; Point(p28) = { -vane_thickness, -ly/2, zmid, cp_vane};
// (COLD) POINTS
p29 = newp; Point(p29) = { +vane_thickness,+ymid, vane_start, cp_vane};
p30 = newp; Point(p30) = { +vane_thickness,+ly/2, zmid, cp_vane};
p31 = newp; Point(p31) = { +vane_thickness,+ymid, lz, cp_vane};
p32 = newp; Point(p32) = { +vane_thickness,-ly/2, zmid, cp_vane};

//Connecting Lines
// Right vane
// HOT
l01 = newl; Line(l01) = {p01, p02};
l02 = newl; Line(l02) = {p02, p03};
l03 = newl; Line(l03) = {p03, p04};
l04 = newl; Line(l04) = {p04, p01};
// COLD
l05 = newl; Line(l05) = {p05, p06};
l06 = newl; Line(l06) = {p06, p07};
l07 = newl; Line(l07) = {p07, p08};
l08 = newl; Line(l08) = {p08, p05};
// DOWN VANE
// HOT
l09 = newl; Line(l09) = {p09, p10};
l10 = newl; Line(l10) = {p10, p11};
l11 = newl; Line(l11) = {p11, p12};
l12 = newl; Line(l12) = {p12, p09};
// COLD
l13 = newl; Line(l13) = {p13, p14};
l14 = newl; Line(l14) = {p14, p15};
l15 = newl; Line(l15) = {p15, p16};
l16 = newl; Line(l16) = {p16, p13};
// LEFT vane
// HOT
l17 = newl; Line(l17) = {p17, p18};
l18 = newl; Line(l18) = {p18, p19};
l19 = newl; Line(l19) = {p19, p20};
l20 = newl; Line(l20) = {p20, p17};
// COLD
l21 = newl; Line(l21) = {p21, p22};
l22 = newl; Line(l22) = {p22, p23};
l23 = newl; Line(l23) = {p23, p24};
l24 = newl; Line(l24) = {p24, p21};
// UP
// HOT
l25 = newl; Line(l25) = {p25, p26};
l26 = newl; Line(l26) = {p26, p27};
l27 = newl; Line(l27) = {p27, p28};
l28 = newl; Line(l28) = {p28, p25};
// COLD
l29 = newl; Line(l29) = {p29, p30};
l30 = newl; Line(l30) = {p30, p31};
l31 = newl; Line(l31) = {p31, p32};
l32 = newl; Line(l32) = {p32, p29};
// RIGHT HOT -> COLD
l33 = newl; Line(l33) = {p01, p05};
l34 = newl; Line(l34) = {p02, p06};
l35 = newl; Line(l35) = {p03, p07};
l36 = newl; Line(l36) = {p04, p08};
// DOWN HOT -> COLD
l37 = newl; Line(l37) = {p09, p13};
l38 = newl; Line(l38) = {p10, p14};
l39 = newl; Line(l39) = {p11, p15};
l40 = newl; Line(l40) = {p12, p16};
// LEFT HOT -> COLD
l41 = newl; Line(l41) = {p17, p21};
l42 = newl; Line(l42) = {p18, p22};
l43 = newl; Line(l43) = {p19, p23};
l44 = newl; Line(l44) = {p20, p24};
// UP HOT -> COLD
l45 = newl; Line(l45) = {p25, p29};
l46 = newl; Line(l46) = {p26, p30};
l47 = newl; Line(l47) = {p27, p31};
l48 = newl; Line(l48) = {p28, p32};

//CURVES
// RIGHT
rh_cl = newcl; Curve Loop(rh_cl) = {l01, l02, l03, l04};
rc_cl = newcl; Curve Loop(rc_cl) = {l05, l06, l07, l08};
rt1_cl = newcl; Curve Loop(rt1_cl) = {+l01,+l34,-l05,-l33};
rt2_cl = newcl; Curve Loop(rt2_cl) = {+l02,+l35,-l06,-l34};
rb1_cl = newcl; Curve Loop(rb1_cl) = {-l04,+l36,+l08,-l33};
rb2_cl = newcl; Curve Loop(rb2_cl) = {-l03,+l35,+l07,-l36};
// DOWN
dh_cl = newcl; Curve Loop(dh_cl) = {l09, l10, l11, l12};
dc_cl = newcl; Curve Loop(dc_cl) = {l13, l14, l15, l16};
dt1_cl = newcl; Curve Loop(dt1_cl) = {+l09,+l38,-l13,-l37};
dt2_cl = newcl; Curve Loop(dt2_cl) = {+l10,+l39,-l14,-l38};
db1_cl = newcl; Curve Loop(db1_cl) = {-l12,+l40,+l16,-l37};
db2_cl = newcl; Curve Loop(db2_cl) = {-l11,+l39,+l15,-l40};
// LEFT
lh_cl = newcl; Curve Loop(lh_cl) = {l17, l18, l19, l20};
lc_cl = newcl; Curve Loop(lc_cl) = {l21, l22, l23, l24};
lt1_cl = newcl; Curve Loop(lt1_cl) = {+l17,+l42,-l21,-l41};
lt2_cl = newcl; Curve Loop(lt2_cl) = {+l18,+l43,-l22,-l42};
lb1_cl = newcl; Curve Loop(lb1_cl) = {-l20,+l44,+l24,-l41};
lb2_cl = newcl; Curve Loop(lb2_cl) = {-l19,+l43,+l23,-l44};
// UP
uh_cl = newcl; Curve Loop(uh_cl) = {l25, l26, l27, l28};
uc_cl = newcl; Curve Loop(uc_cl) = {l29, l30, l31, l32};
ut_cl = newcl; Curve Loop(ut_cl) = {+l25, +l26, +l47, -l30, -l29, -l45};
ub_cl = newcl; Curve Loop(ub_cl) = {-l28, -l27, +l47, +l31, +l32, -l45};
ut1_cl = newcl; Curve Loop(ut1_cl) = {+l25,+l46,-l29,-l45};
ut2_cl = newcl; Curve Loop(ut2_cl) = {+l26,+l47,-l30,-l46};
ub1_cl = newcl; Curve Loop(ub1_cl) = {-l28,+l48,+l32,-l45};
ub2_cl = newcl; Curve Loop(ub2_cl) = {-l27,+l47,+l31,-l48};
// Surfaces
// RIGHT VANE
rh_s = news; Surface(rh_s) = {rh_cl};
rc_s = news; Surface(rc_s) = {rc_cl};
rt1_s = news; Surface(rt1_s) = {rt1_cl};
rb1_s = news; Surface(rb1_s) = {rb1_cl};
rt2_s = news; Surface(rt2_s) = {rt2_cl};
rb2_s = news; Surface(rb2_s) = {rb2_cl};
// DOWN VANE
dh_s = news; Surface(dh_s) = {dh_cl};
dc_s = news; Surface(dc_s) = {dc_cl};
dt1_s = news; Surface(dt1_s) = {dt1_cl};
db1_s = news; Surface(db1_s) = {db1_cl};
dt2_s = news; Surface(dt2_s) = {dt2_cl};
db2_s = news; Surface(db2_s) = {db2_cl};
// LEFT VANE
lh_s = news; Surface(lh_s) = {lh_cl};
lc_s = news; Surface(lc_s) = {lc_cl};
lt1_s = news; Surface(lt1_s) = {lt1_cl};
lb1_s = news; Surface(lb1_s) = {lb1_cl};
lt2_s = news; Surface(lt2_s) = {lt2_cl};
lb2_s = news; Surface(lb2_s) = {lb2_cl};
// UP VANE
uh_s = news; Surface(uh_s) = {uh_cl};
uc_s = news; Surface(uc_s) = {uc_cl};
ut1_s = news; Surface(ut1_s) = {ut1_cl};
ub1_s = news; Surface(ub1_s) = {ub1_cl};
ut2_s = news; Surface(ut2_s) = {ut2_cl};
ub2_s = news; Surface(ub2_s) = {ub2_cl};
//
// Physical Surfaces
// RIGHT VANE
Physical Surface("RightHot",188) = {rh_s};
Physical Surface("RightCold",183) = {rc_s};
Physical Surface("RightTransition",1820) = {rt1_s, rt2_s,rb1_s, rb2_s};
//Physical Surface("RightTransitionTop",182020) = {rt1_s, rt2_s};
//Physical Surface("RightTransitionBottom",18202) = {rb1_s, rb2_s};
// DOWN VANE
Physical Surface("DownHot", 48) = {dh_s};
Physical Surface("DownCold",43) = {dc_s};
Physical Surface("DownTransition",420) = {dt1_s, dt2_s, db1_s, db2_s};
//Physical Surface("DownTransitionTop",42020) = {dt1_s, dt2_s};
//Physical Surface("DownTransitionBottom",4202) = {db1_s, db2_s};
// LEFT VANE
Physical Surface("LeftHot",128) = {lh_s};
Physical Surface("LeftCold",123) = {lc_s};
Physical Surface("LeftTransition",1220) = {lt1_s, lt2_s, lb1_s, lb2_s};
//Physical Surface("LeftTransitionTop",122020) = {lt1_s, lt2_s};
//Physical Surface("LeftTransitionBottom",12202) = {lb1_s, lb2_s};
// UP VANE
Physical Surface("UpHot", 218) = {uh_s};
Physical Surface("UpCold",213) = {uc_s};
Physical Surface("UpTransition",2120) = {ut1_s, ut2_s, ub1_s, ub2_s};
//Physical Surface("UpTransitionTop",212020) = {ut1_s, ut2_s};
//Physical Surface("UpTransitionBottom",21202) = {ub1_s, ub2_s};
//+
r_sl = newsl; Surface Loop(r_sl) = {rh_s, rc_s, rt1_s, rt2_s, rb1_s, rb2_s};
d_sl = newsl; Surface Loop(d_sl) = {dh_s, dc_s, dt1_s, dt2_s, db1_s, db2_s};
l_sl = newsl; Surface Loop(l_sl) = {lh_s, lc_s, lt1_s, lt2_s, lb1_s, lb2_s};
u_sl = newsl; Surface Loop(u_sl) = {uh_s, uc_s, ut1_s, ut2_s, ub1_s, ub2_s};
v_r =newv; Volume(v_r) = {r_sl};
v_d =newv; Volume(v_d) = {d_sl};
v_l =newv; Volume(v_l) = {l_sl};
v_u =newv; Volume(v_u) = {u_sl};
//+
glass = newv; Sphere(glass) = {0, 0, 0, radius};
//
Physical Surface("Wall", 23) = {99};
//
 Physical Volume ("VOLUME",  8888) = { BooleanDifference{ Volume{glass}; Delete;
 }{ Volume{v_r}; Volume{v_d}; Volume{v_l}; Volume{v_u}; }
 };
 
