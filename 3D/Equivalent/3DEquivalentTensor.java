/// Remove old geometry if it exists
String[] geomTags = model.geom().tags();
for (String tag : geomTags) {
  if (tag.equals("geom1")) {
    model.geom().remove(tag);
    break;
  }
}



/// Create new 3D geometry
model.geom().create("geom1", 3);
model.geom("geom1").lengthUnit("m");



/// Define parameters
String stdn = "1"; // Number of the new study
//lattice sizing
int nucx = 4; //number of unit cells in x
int nucy = 4; //number of unit cells in y
int nucz = 4; //number of unit cells in z
double[] block_sides_outer = {0.08, 0.08, 0.08}; //size of unit cells (x,y,z)
double avgside = (block_sides_outer[0]+block_sides_outer[1]+block_sides_outer[2])/3;
//unit cell settings
double radius = 0.004; //radius of poles (for 0.08, 0.0012 is good min)
double rns = radius/avgside;
double wall_thickness = radius; //thickenss of walls
//inner block stuff
double[] theta = {0, 0};

//fix - actually get the value or two values
double kd = 238; // Real isotropic conductivity value of aluminium 3003
double E = -43520*rns*rns*rns*rns+11610*rns*rns*rns-1166*rns*rns+57.02*rns-0.6073;
if (rns > 0.0875) {
  E = -35.5*rns*rns*rns*rns+19.8*rns*rns*rns+4.196*rns*rns+1.528*rns+0.5088;
  if (rns > 0.1375) {
    E = -602.6*rns*rns*rns*rns+532.4*rns*rns*rns-178*rns*rns+27.62*rns-0.7644;
    if (rns > 0.2875) {
      E = 0.02667*rns+0.9913;
      if (rns > 0.325) {
        E = 1;
      }
    }
  }
}



///Input block-range for theta[0] is 0 to 90, range for theta[1] is anything(normally-180 to 180)
double[][][][] block = {
  {{{10, -135}, {10, -100}, {10, -80}, {10, -45}},
   {{10, -170}, {45, -135}, {45, -45}, {10, -10}},
   {{10, 170}, {45, 135}, {45, 45}, {10, 10}},
   {{10, 135}, {10, 100}, {10, 80}, {10, 45}}},
  
  {{{35, -135}, {35, -100}, {35, -80}, {35, -45}},
   {{35, -170}, {85, -135}, {85, -45}, {35, -10}},
   {{35, 170}, {85, 135}, {85, 45}, {35, 10}},
   {{35, 135}, {35, 100}, {35, 80}, {35, 45}}},
  
  {{{45, -135}, {45, -100}, {45, -80}, {45, -45}},
   {{45, -170}, {10, -135}, {10, -45}, {45, -10}},
   {{45, 170}, {10, 135}, {10, 45}, {45, 10}},
   {{45, 135}, {45, 100}, {45, 80}, {45, 45}}},
  
  {{{75, -135}, {75, -100}, {75, -80}, {75, -45}},
   {{75, -170}, {30, -135}, {30, -45}, {75, -10}},
   {{75, 170}, {30, 135}, {30, 45}, {75, 10}},
   {{75, 135}, {75, 100}, {75, 80}, {75, 45}}}
};



/// Loop through block
for (int i = 0; i < nucz; i++) {
  for (int j = 0; j < nucy; j++) {
    for (int k = 0; k < nucx; k++) {
      theta[0] = block[i][j][k][0];
      theta[1] = block[i][j][k][1];
      String tag = String.format("bl_%d_%d_%d", i, j, k);
      double[] ref_fr = {k*block_sides_outer[0], j*block_sides_outer[1], i*block_sides_outer[2]}; // Reference frame
      
      if (theta[0] == 360) {
        continue;
      }
      
      /// 4 Build outside block
      String outerTag = tag+"_outer";
      model.geom("geom1").feature().create(outerTag, "Block");
      model.geom("geom1").feature(outerTag).set("base", "corner");
      model.geom("geom1").feature(outerTag).set("size", block_sides_outer);
      model.geom("geom1").feature(outerTag).set("pos", ref_fr);
      
      
      
      // Calculate effective conductivity matrix
      double c1 = Math.cos(theta[0]*Math.PI/180);
      double c2 = Math.cos(theta[1]*Math.PI/180);
      double s1 = Math.sin(theta[0]*Math.PI/180);
      double s2 = Math.sin(theta[1]*Math.PI/180);
      double kxx = s1*s1*c2*c2+E*c1*c1*c2*c2+E*s2*s2;
      double kyy = s1*s1*s2*s2+E*c1*c1*s2*s2+E*c2*c2;
      double kzz = c1*c1+E*s1*s1;
      double kyx = (s1*s1+E*c1*c1-E)*s2*c2;
      double kzx = (1-E)*s1*c1*c2;
      double kzy = (1-E)*s1*c1*s2;
      //symmetry on
      double kxy = kyx;
      double kxz = kzx;
      double kyz = kzy;
      
      String sm = "sm"+i+0000+j+000+k;
      model.selection().create(sm, "Box");
      model.selection(sm).geom("geom1", 3);
      model.selection(sm).set("xmin", ref_fr[0]+radius*0.1);
      model.selection(sm).set("xmax", ref_fr[0]+radius*0.2);
      model.selection(sm).set("ymin", ref_fr[1]+radius*0.1);
      model.selection(sm).set("ymax", ref_fr[1]+radius*0.2);
      model.selection(sm).set("zmin", ref_fr[2]+radius*0.1);
      model.selection(sm).set("zmax", ref_fr[2]+radius*0.2);
      
      String mat = "mat"+i+0000+j+0000+k;
      model.material().create(mat, "Common", "mod1");
      model.material(mat).label(mat);
      with(model.material(mat).propertyGroup("def"));
        set("thermalconductivity", new String[]{
          kd*kxx+"[W/(m*K)]", kd*kxy+"[W/(m*K)]", kd*kxz+"[W/(m*K)]",
          kd*kyx+"[W/(m*K)]", kd*kyy+"[W/(m*K)]", kd*kyz+"[W/(m*K)]",
          kd*kzx+"[W/(m*K)]", kd*kzy+"[W/(m*K)]", kd*kzz+"[W/(m*K)]"
        });
        set("density", new String[]{"2700[kg/m^3]"});
        set("heatcapacity", new String[]{"900[J/(kg*K)]"});
      endwith();
      model.material(mat).selection().named(sm);
    }
  }
}



/// Generate Geometry
model.geom("geom1").run();



/// Calculations for physics
//hot
for (int i = 0; i < nucx; i++) {
  for (int j = 0; j < nucy; j++) {
    String sel_hotsurface = "hsf"+i+j*0.001;
    double xli = i*block_sides_outer[0]+radius;
    double yli = j*block_sides_outer[1]+radius;
    double zli = 0;
    model.selection().create(sel_hotsurface, "Box");
    model.selection(sel_hotsurface).geom("geom1", 2);
    model.selection(sel_hotsurface).set("xmin", xli);
    model.selection(sel_hotsurface).set("xmax", xli);
    model.selection(sel_hotsurface).set("ymin", yli);
    model.selection(sel_hotsurface).set("ymax", yli);
    model.selection(sel_hotsurface).set("zmin", zli);
    model.selection(sel_hotsurface).set("zmax", zli);
  }
}
//cold
for (int i = 0; i < nucx; i++) {
  for (int j = 0; j < nucy; j++) {
    String sel_coldsurface = "csf"+i+j*0.001;
    double xli = i*block_sides_outer[0]+radius;
    double yli = j*block_sides_outer[1]+radius;
    double zli = nucz*block_sides_outer[2];
    model.selection().create(sel_coldsurface, "Box");
    model.selection(sel_coldsurface).geom("geom1", 2);
    model.selection(sel_coldsurface).set("xmin", xli);
    model.selection(sel_coldsurface).set("xmax", xli+0.1*radius);
    model.selection(sel_coldsurface).set("ymin", yli);
    model.selection(sel_coldsurface).set("ymax", yli+0.1*radius);
    model.selection(sel_coldsurface).set("zmin", zli-0.1*radius);
    model.selection(sel_coldsurface).set("zmax", zli+0.1*radius);
  }
}




/// Add physics
model.physics().create("ht", "HeatTransfer", "geom1");
//hot
for (int i = 0; i < nucx; i++) {
  for (int j = 0; j < nucy; j++) {
    String sel_hotsurface = "hsf"+i+j*0.001;
    String temp = "temph"+i+j*0.001;
    model.physics("ht").create(temp, "TemperatureBoundary", 2);
    model.physics("ht").feature(temp).selection().named(sel_hotsurface);
    model.physics("ht").feature(temp).set("T0", "400[K]");
  }
}
//cold
for (int i = 0; i < nucx; i++) {
  for (int j = 0; j < nucy; j++) {
    String sel_coldsurface = "csf"+i+j*0.001;
    String temp = "tempc"+i+j*0.001;
    model.physics("ht").create(temp, "TemperatureBoundary", 2);
    model.physics("ht").feature(temp).selection().named(sel_coldsurface);
    model.physics("ht").feature(temp).set("T0", "300[K]");
  }
}



/// Add mesh
model.mesh().create("mesh1", "geom1");
model.mesh("mesh1").automatic(true); // Keep on Physics-controlled mesh
model.mesh("mesh1").autoMeshSize(6); // [1,9] [fine,coarse]
model.mesh("mesh1").run();



/// Add study
model.study().create(stdn);
model.study(stdn).create("stat", "Stationary");
model.study(stdn).run();
