//+
SetFactory("OpenCASCADE");
//Control points
Mesh.CharacteristicLengthMax = 0.1;
// inner control points
cp_vane = 0.01;
//
Mesh.RandomFactor = 5.e-10;
// Choose 2D Algorithm
// 1 : MeshAdapt
// 2 : Automatic
// 5 : Delaunay
// 6 : Frontal-Delaunay
Mesh.Algorithm = 1;
//Choose 3D Algorithm
// 1 : Delaunay
Mesh.Algorithm3D =1;

vane_length = 1; 
xmid = 0.75 * vane_length;
ymid = 0;
zmid = 0.75 * vane_length;
//radius = 2 * vane_length;
radius = 2.5;
vane_thickness = 0.05; //total vane thickness 0.05*2
vane_start = 0.25; //xmid - (vane_length/2);
lx = vane_start + vane_length;
ly = vane_length;
lz = vane_start + vane_length;
//
// Vane RIGHT
// Define points - TOP
p_rt1 = newp;  Point(p_rt1) = {+vane_start, +ly/2, +vane_thickness, cp_vane};
p_rt2 = newp;  Point(p_rt2) = {+lx, +ly/2, +vane_thickness, cp_vane};
p_rt3 = newp;  Point(p_rt3) = {+lx, +ly/2, -vane_thickness, cp_vane};
p_rt4 = newp;  Point(p_rt4) = {+vane_start, +ly/2, -vane_thickness, cp_vane};
// Define points - Bottom
p_rb1 = newp;  Point(p_rb1) = {+vane_start, -ly/2, +vane_thickness, cp_vane};
p_rb2 = newp;  Point(p_rb2) = {+lx, -ly/2, +vane_thickness, cp_vane};
p_rb3 = newp;  Point(p_rb3) = {+lx, -ly/2, -vane_thickness, cp_vane};
p_rb4 = newp;  Point(p_rb4) = {+vane_start, -ly/2, -vane_thickness, cp_vane};
// Connect Lines - Top
l_rt1 = newl; Line(l_rt1) = {p_rt1, p_rt2};
l_rt2 = newl; Line(l_rt2) = {p_rt2, p_rt3};
l_rt3 = newl; Line(l_rt3) = {p_rt3, p_rt4};
l_rt4 = newl; Line(l_rt4) = {p_rt4, p_rt1};
// Connect Lines - Bottom 
l_rb1 = newl; Line(l_rb1) = {p_rb1, p_rb2};
l_rb2 = newl; Line(l_rb2) = {p_rb2, p_rb3};
l_rb3 = newl; Line(l_rb3) = {p_rb3, p_rb4};
l_rb4 = newl; Line(l_rb4) = {p_rb4, p_rb1};
// Connect Lines - Top and bottom
l_rtb1 = newl; Line(l_rtb1) = {p_rt1, p_rb1};
l_rtb2 = newl; Line(l_rtb2) = {p_rt2, p_rb2};
l_rtb3 = newl; Line(l_rtb3) = {p_rt3, p_rb3};
l_rtb4 = newl; Line(l_rtb4) = {p_rt4, p_rb4};
// Vane Left
// Define points - TOP
p_lt1 = newp;  Point(p_lt1) = {-vane_start, +ly/2, +vane_thickness, cp_vane};
p_lt2 = newp;  Point(p_lt2) = {-lx, +ly/2, +vane_thickness, cp_vane};
p_lt3 = newp;  Point(p_lt3) = {-lx, +ly/2, -vane_thickness, cp_vane};
p_lt4 = newp;  Point(p_lt4) = {-vane_start, +ly/2, -vane_thickness, cp_vane};
// Define points - Bottom
p_lb1 = newp;  Point(p_lb1) = {-vane_start, -ly/2, +vane_thickness, cp_vane};
p_lb2 = newp;  Point(p_lb2) = {-lx, -ly/2, +vane_thickness, cp_vane};
p_lb3 = newp;  Point(p_lb3) = {-lx, -ly/2, -vane_thickness, cp_vane};
p_lb4 = newp;  Point(p_lb4) = {-vane_start, -ly/2, -vane_thickness, cp_vane};
// Connect Lines - Top
l_lt1 = newl; Line(l_lt1) = {p_lt1, p_lt2};
l_lt2 = newl; Line(l_lt2) = {p_lt2, p_lt3};
l_lt3 = newl; Line(l_lt3) = {p_lt3, p_lt4};
l_lt4 = newl; Line(l_lt4) = {p_lt4, p_lt1};
// Connect Lines - Bottom 
l_lb1 = newl; Line(l_lb1) = {p_lb1, p_lb2};
l_lb2 = newl; Line(l_lb2) = {p_lb2, p_lb3};
l_lb3 = newl; Line(l_lb3) = {p_lb3, p_lb4};
l_lb4 = newl; Line(l_lb4) = {p_lb4, p_lb1};
// Connect Lines - Top and bottom
l_ltb1 = newl; Line(l_ltb1) = {p_lt1, p_lb1};
l_ltb2 = newl; Line(l_ltb2) = {p_lt2, p_lb2};
l_ltb3 = newl; Line(l_ltb3) = {p_lt3, p_lb3};
l_ltb4 = newl; Line(l_ltb4) = {p_lt4, p_lb4};
// Vane DOWN
// Define points - TOP
p_dt1 = newp;  Point(p_dt1) = {+vane_thickness, +ly/2, -vane_start, cp_vane};
p_dt2 = newp;  Point(p_dt2) = {+vane_thickness, +ly/2, -lz, cp_vane};
p_dt3 = newp;  Point(p_dt3) = {-vane_thickness, +ly/2, -lz, cp_vane};
p_dt4 = newp;  Point(p_dt4) = {-vane_thickness, +ly/2, -vane_start, cp_vane};
// Define points - Bottom
p_db1 = newp;  Point(p_db1) = {+vane_thickness, -ly/2, -vane_start, cp_vane};
p_db2 = newp;  Point(p_db2) = {+vane_thickness, -ly/2, -lz, cp_vane};
p_db3 = newp;  Point(p_db3) = {-vane_thickness, -ly/2, -lz, cp_vane};
p_db4 = newp;  Point(p_db4) = {-vane_thickness, -ly/2, -vane_start, cp_vane};
// Connect Lines - Top
l_dt1 = newl; Line(l_dt1) = {p_dt1, p_dt2};
l_dt2 = newl; Line(l_dt2) = {p_dt2, p_dt3};
l_dt3 = newl; Line(l_dt3) = {p_dt3, p_dt4};
l_dt4 = newl; Line(l_dt4) = {p_dt4, p_dt1};
// Connect Lines - Bottom 
l_db1 = newl; Line(l_db1) = {p_db1, p_db2};
l_db2 = newl; Line(l_db2) = {p_db2, p_db3};
l_db3 = newl; Line(l_db3) = {p_db3, p_db4};
l_db4 = newl; Line(l_db4) = {p_db4, p_db1};
// Connect Lines - Top and bottom
l_dtb1 = newl; Line(l_dtb1) = {p_dt1, p_db1};
l_dtb2 = newl; Line(l_dtb2) = {p_dt2, p_db2};
l_dtb3 = newl; Line(l_dtb3) = {p_dt3, p_db3};
l_dtb4 = newl; Line(l_dtb4) = {p_dt4, p_db4};
// Vane Up 
// Define points - TOP
p_ut1 = newp;  Point(p_ut1) = {+vane_thickness, +ly/2, +vane_start, cp_vane};
p_ut2 = newp;  Point(p_ut2) = {+vane_thickness, +ly/2, +lz, cp_vane};
p_ut3 = newp;  Point(p_ut3) = {-vane_thickness, +ly/2, +lz, cp_vane};
p_ut4 = newp;  Point(p_ut4) = {-vane_thickness, +ly/2, +vane_start, cp_vane};
// Define points - Bottom
p_ub1 = newp;  Point(p_ub1) = {+vane_thickness, -ly/2, +vane_start, cp_vane};
p_ub2 = newp;  Point(p_ub2) = {+vane_thickness, -ly/2, +lz, cp_vane};
p_ub3 = newp;  Point(p_ub3) = {-vane_thickness, -ly/2, +lz, cp_vane};
p_ub4 = newp;  Point(p_ub4) = {-vane_thickness, -ly/2, +vane_start, cp_vane};
// Connect Lines - Top
l_ut1 = newl; Line(l_ut1) = {p_ut1, p_ut2};
l_ut2 = newl; Line(l_ut2) = {p_ut2, p_ut3};
l_ut3 = newl; Line(l_ut3) = {p_ut3, p_ut4};
l_ut4 = newl; Line(l_ut4) = {p_ut4, p_ut1};
// Connect Lines - Bottom 
l_ub1 = newl; Line(l_ub1) = {p_ub1, p_ub2};
l_ub2 = newl; Line(l_ub2) = {p_ub2, p_ub3};
l_ub3 = newl; Line(l_ub3) = {p_ub3, p_ub4};
l_ub4 = newl; Line(l_ub4) = {p_ub4, p_ub1};
// Connect Lines - Top and bottom
l_utb1 = newl; Line(l_utb1) = {p_ut1, p_ub1};
l_utb2 = newl; Line(l_utb2) = {p_ut2, p_ub2};
l_utb3 = newl; Line(l_utb3) = {p_ut3, p_ub3};
l_utb4 = newl; Line(l_utb4) = {p_ut4, p_ub4};

// Curves Loops TOP of each vane 
//right, down, left, up
cl_r_top = newcl;Curve Loop(cl_r_top) = {l_rt1, l_rt2, -l_rt3, l_rt4};
cl_d_top = newcl;Curve Loop(cl_d_top) = {l_dt1, l_dt2, -l_dt3, l_dt4};
cl_l_top = newcl;Curve Loop(cl_l_top) = {l_lt1, l_lt2, -l_lt3, -l_lt4};
cl_u_top = newcl;Curve Loop(cl_u_top) = {l_ut1, l_ut2, -l_ut3, -l_ut4};
// Curves Loops BOTTOM of each vane 
//right, down, left, up
cl_r_bot = newcl;Curve Loop(cl_r_bot) = {l_rb1, l_rb2, -l_rb3, l_rb4};
cl_d_bot = newcl;Curve Loop(cl_d_bot) = {l_db1, l_db2, -l_db3, l_db4};
cl_l_bot = newcl;Curve Loop(cl_l_bot) = {l_lb1, l_lb2, -l_lb3, -l_lb4};
cl_u_bot = newcl;Curve Loop(cl_u_bot) = {l_ub1, l_ub2, -l_ub3, -l_ub4};
// Curves Loops Hot-side & COLD-side of each vane 
//right, down, left, up
cl_r_hot = newcl;Curve Loop(cl_r_hot) = {l_rt1, l_rtb2, -l_rb1, -l_rtb1};
cl_r_col = newcl;Curve Loop(cl_r_col) = {-l_rt3, l_rtb3, l_rb3, -l_rtb4};

cl_d_hot = newcl;Curve Loop(cl_d_hot) = {l_dt1, l_dtb2, -l_db1, -l_dtb1};
cl_d_col = newcl;Curve Loop(cl_d_col) = {-l_dt3, l_dtb3, l_db3, -l_dtb4};

cl_l_col = newcl;Curve Loop(cl_l_col) = {l_lt1, l_ltb2, -l_lb1, -l_ltb1};
cl_l_hot = newcl;Curve Loop(cl_l_hot) = {-l_lt3, l_ltb3, l_lb3, -l_ltb4};

cl_u_col = newcl;Curve Loop(cl_u_col) = {l_ut1, l_utb2, -l_ub1, -l_utb1};
cl_u_hot = newcl;Curve Loop(cl_u_hot) = {-l_ut3, l_utb3, l_ub3, -l_utb4};
//right, down, left, up
// Curves Loops outer-Side (is) of each vane 
//right, down, left, up
cl_r_os = newcl;Curve Loop(cl_r_os) = {l_rt2, l_rtb3, -l_rb2, -l_rtb2};
cl_d_os = newcl;Curve Loop(cl_d_os) = {l_dt2, l_dtb3, -l_db2, -l_dtb2};
cl_l_os = newcl;Curve Loop(cl_l_os) = {l_lt2, l_ltb3, -l_lb2, -l_ltb2};
cl_u_os = newcl;Curve Loop(cl_u_os) = {l_ut2, l_utb3, -l_ub2, -l_utb2};
// Curves Loops inner-Side (os) of each vane 
//right, down, left, up
cl_r_is = newcl;Curve Loop(cl_r_is) = {l_rt4, l_rtb1, -l_rb4, -l_rtb4};
cl_d_is = newcl;Curve Loop(cl_d_is) = {l_dt4, l_dtb1, -l_db4, -l_dtb4};
cl_l_is = newcl;Curve Loop(cl_l_is) = {l_lt4, l_ltb1, -l_lb4, -l_ltb4};
cl_u_is = newcl;Curve Loop(cl_u_is) = {l_ut4, l_utb1, -l_ub4, -l_utb4};
// Construct Suface of each vane
// Top surface - right, down, left, right
s_r_top = news; Surface(s_r_top)  = {cl_r_top};
s_d_top = news; Surface(s_d_top)  = {cl_d_top};
s_l_top = news; Surface(s_l_top)  = {cl_l_top};
s_u_top = news; Surface(s_u_top)  = {cl_u_top};
// Bottom surface - right, down, left, right
s_r_bot = news; Surface(s_r_bot)  = {cl_r_bot};
s_d_bot = news; Surface(s_d_bot)  = {cl_d_bot};
s_l_bot = news; Surface(s_l_bot)  = {cl_l_bot};
s_u_bot = news; Surface(s_u_bot)  = {cl_u_bot};
// Inner-side surface - right, down, left, right
s_r_is = news; Surface(s_r_is)  = {cl_r_is};
s_d_is = news; Surface(s_d_is)  = {cl_d_is};
s_l_is = news; Surface(s_l_is)  = {cl_l_is};
s_u_is = news; Surface(s_u_is)  = {cl_u_is};
// Outer-side surface - right, down, left, right
s_r_os = news; Surface(s_r_os)  = {cl_r_os};
s_d_os = news; Surface(s_d_os)  = {cl_d_os};
s_l_os = news; Surface(s_l_os)  = {cl_l_os};
s_u_os = news; Surface(s_u_os)  = {cl_u_os};
// Hot surface - right, down, left, right
s_r_hot = news; Surface(s_r_hot)  = {cl_r_hot};
s_d_hot = news; Surface(s_d_hot)  = {cl_d_hot};
s_l_hot = news; Surface(s_l_hot)  = {cl_l_hot};
s_u_hot = news; Surface(s_u_hot)  = {cl_u_hot};
// Cold surface - right, down, left, right
s_r_col = news; Surface(s_r_col)  = {cl_r_col};
s_d_col = news; Surface(s_d_col)  = {cl_d_col};
s_l_col = news; Surface(s_l_col)  = {cl_l_col};
s_u_col = news; Surface(s_u_col)  = {cl_u_col};

Physical Surface("RIGHTCOLD", 183) = {s_r_col};
Physical Surface("DOWNCOLD", 43) = {s_d_col};
Physical Surface("LEFTCOLD", 123) = {s_l_col};
Physical Surface("UPCOLD", 213) = {s_u_col};
Physical Surface("RIGHTHOT", 188) = {s_r_hot};
Physical Surface("DOWNHOT", 48) = {s_d_hot};
Physical Surface("LEFTHOT", 128) = {s_l_hot};
Physical Surface("UPHOT", 218) = {s_u_hot};
Physical Surface("RightTransition", 1820) = {s_r_top, s_r_bot, s_r_is, s_r_os};
//Physical Surface("RightTransitionTop", 182020) = {s_r_top};
//Physical Surface("RightTransitionBot", 18202) = {s_r_bot};
//Physical Surface("RightTransitionInner", 18209) = {s_r_is};
//Physical Surface("RightTransitionOuter", 182015) = {s_r_os};
//
Physical Surface("DownTransition", 420) = {s_d_top, s_d_bot, s_d_is, s_d_os};
//Physical Surface("DownTransitionTop", 42020) = {s_d_top};
//Physical Surface("DownTransitionBot", 4202) = {s_d_bot};
//Physical Surface("DownTransitionInner", 4209) = {s_d_is};
//Physical Surface("DownTransitionOuter", 42015) = {s_d_os};
//
Physical Surface("LeftTransition", 1220) = {s_l_top, s_l_bot, s_l_is, s_l_os};
//Physical Surface("LeftTransitionTop", 122020) = {s_l_top};
//Physical Surface("LeftTransitionBot", 12202) = {s_l_bot};
//Physical Surface("LeftTransitionInner", 12209) = {s_l_is};
//Physical Surface("LeftTransitionOuter", 122015) = {s_l_os};
//
Physical Surface("UpTransition", 2120) = {s_u_top, s_u_bot, s_u_is, s_u_os};
//Physical Surface("UpTransitionTop", 212020) = {s_u_top};
//Physical Surface("UpTransitionBot", 21202) = {s_u_bot};
//Physical Surface("UpTransitionInner", 21209) = {s_u_is};
//Physical Surface("UpTransitionOut", 212015) = {s_u_os};
// Surface Loop
sl_r  = newsl; Surface Loop(sl_r) = { s_r_top, s_r_os, s_r_bot, s_r_is, s_r_col, s_r_hot };
sl_d  = newsl; Surface Loop(sl_d) = { s_d_top, s_d_os, s_d_bot, s_d_is, s_d_col, s_d_hot };
sl_l  = newsl; Surface Loop(sl_l) = { s_l_top, s_l_os, s_l_bot, s_l_is, s_l_col, s_l_hot };
sl_u  = newsl; Surface Loop(sl_u) = { s_u_top, s_u_os, s_u_bot, s_u_is, s_u_col, s_u_hot };
// Vane
v_r = newv; Volume(v_r) = {sl_r};
//
v_d = newv; Volume(v_d) = {sl_d};
//
v_l = newv; Volume(v_l) = {sl_l};
//
v_u = newv; Volume(v_u) = {sl_u};
//+
glass = newv; Sphere(glass) = {0, 0, 0, radius};
// Get the id of the newly created volume
//Printf("%g", glass[0]);
Physical Surface("Wall", 23) = {97};
//+
Physical Volume ("VOLUME",8888) = { BooleanDifference{ Volume{glass}; Delete; }{ Volume{v_r}; Volume{v_d}; Volume{v_l}; Volume{v_u};} };
