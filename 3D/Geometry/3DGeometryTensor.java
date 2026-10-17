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
//unit cell settings
double radius = 0.004; //radius of poles
double wall_thickness = radius; //thickenss of walls
//inside block
double[] block_sides_inner = {block_sides_outer[0]-2*wall_thickness, block_sides_outer[1]-2*wall_thickness, block_sides_outer[2]-2*wall_thickness}; //size of cutout block
//inner block stuff
double[] theta = {0, 0};



/// Input block - range for theta[0] is 0 to 90, range for theta[1] is anything (normally -180 to 180)
double[][][][] block = {
  {{{10, -135}, {10, -115}, {10, -100}, {10, -80}, {10, -65}, {10, -45}},
   {{10, -155}, {25, -135}, {25, -115}, {25, -65}, {25, -45}, {10, -25}},
   {{10, -170}, {25, -155}, {45, -135}, {45, -45}, {25, -25}, {10, -10}},
   {{10, 170}, {25, 155}, {45, 135}, {45, 45}, {25, 25}, {10, 10}},
   {{10, 155}, {25, 135}, {25, 115}, {25, 65}, {25, 45}, {10, 25}},
   {{10, 135}, {10, 115}, {10, 100}, {10, 80}, {10, 65}, {10, 45}}},
  
  {{{20, -135}, {20, -115}, {20, -100}, {20, -80}, {20, -65}, {20, -45}},
   {{20, -155}, {45, -135}, {45, -115}, {45, -65}, {45, -45}, {20, -25}},
   {{20, -170}, {45, -155}, {60, -135}, {60, -45}, {45, -25}, {20, -10}},
   {{20, 170}, {45, 155}, {60, 135}, {60, 45}, {45, 25}, {20, 10}},
   {{20, 155}, {45, 135}, {45, 115}, {45, 65}, {45, 45}, {20, 25}},
   {{20, 135}, {20, 115}, {20, 100}, {20, 80}, {20, 65}, {20, 45}}},
  
  {{{35, -135}, {35, -115}, {35, -100}, {35, -80}, {35, -65}, {35, -45}},
   {{35, -155}, {70, -135}, {70, -115}, {70, -65}, {70, -45}, {35, -25}},
   {{35, -170}, {70, -155}, {85, -135}, {85, -45}, {70, -25}, {35, -10}},
   {{35, 170}, {70, 155}, {85, 135}, {85, 45}, {70, 25}, {35, 10}},
   {{35, 155}, {70, 135}, {70, 115}, {70, 65}, {70, 45}, {35, 25}},
   {{35, 135}, {35, 115}, {35, 100}, {35, 80}, {35, 65}, {35, 45}}},
  
  {{{45, -135}, {45, -115}, {45, -100}, {45, -80}, {45, -65}, {45, -45}},
   {{45, -155}, {20, -135}, {20, -115}, {20, -65}, {20, -45}, {45, -25}},
   {{45, -170}, {20, -155}, {10, -135}, {10, -45}, {20, -25}, {45, -10}},
   {{45, 170}, {20, 155}, {10, 135}, {10, 45}, {20, 25}, {45, 10}},
   {{45, 155}, {20, 135}, {20, 115}, {20, 65}, {20, 45}, {45, 25}},
   {{45, 135}, {45, 115}, {45, 100}, {45, 80}, {45, 65}, {45, 45}}},
  
  {{{60, -135}, {60, -115}, {60, -100}, {60, -80}, {60, -65}, {60, -45}},
   {{60, -155}, {30, -135}, {30, -115}, {30, -65}, {30, -45}, {60, -25}},
   {{60, -170}, {30, -155}, {20, -135}, {20, -45}, {30, -25}, {60, -10}},
   {{60, 170}, {30, 155}, {20, 135}, {20, 45}, {30, 25}, {60, 10}},
   {{60, 155}, {30, 135}, {30, 115}, {30, 65}, {30, 45}, {60, 25}},
   {{60, 135}, {60, 115}, {60, 100}, {60, 80}, {60, 65}, {60, 45}}},
  
  {{{75, -135}, {75, -115}, {75, -100}, {75, -80}, {75, -65}, {75, -45}},
   {{75, -155}, {40, -135}, {40, -115}, {40, -65}, {40, -45}, {75, -25}},
   {{75, -170}, {40, -155}, {30, -135}, {30, -45}, {40, -25}, {75, -10}},
   {{75, 170}, {40, 155}, {30, 135}, {30, 45}, {40, 25}, {75, 10}},
   {{75, 155}, {40, 135}, {40, 115}, {40, 65}, {40, 45}, {75, 25}},
   {{75, 135}, {75, 115}, {75, 100}, {75, 80}, {75, 65}, {75, 45}}}
};



/// Loop through block
for (int i = 0; i < nucz; i++) {
  for (int j = 0; j < nucy; j++) {
    for (int k = 0; k < nucx; k++) {
      theta[0] = block[i][j][k][0];
      theta[1] = block[i][j][k][1];
      //theta[0] = 45;
      //theta[1] = 45;
      String tag = String.format("bl_%d_%d_%d", i, j, k);
      double[] ref_fr = {k*block_sides_outer[0], j*block_sides_outer[1], i*block_sides_outer[2]}; // Reference frame
      
      
      
      /// 1 Build temporary inside block
      double[] block_inner = {ref_fr[0]+wall_thickness, ref_fr[1]+wall_thickness, ref_fr[2]+wall_thickness};
      String tempinnerTag = tag+"_temp_inner";
      model.geom("geom1").feature().create(tempinnerTag, "Block");
      model.geom("geom1").feature(tempinnerTag).set("base", "corner");
      model.geom("geom1").feature(tempinnerTag).set("size", block_sides_inner);
      model.geom("geom1").feature(tempinnerTag).set("pos", block_inner);
      
      /// 2 Build poles
      double[] tilt = {Math.sin(theta[0]*Math.PI/180)*Math.cos(theta[1]*Math.PI/180), Math.sin(theta[0]*Math.PI/180)*Math.sin(theta[1]*Math.PI/180), Math.sqrt(1-Math.sin(theta[0]*Math.PI/180)*Math.sin(theta[0]*Math.PI/180))}; // Calculate unit vector of poles
      double shift = 0.5*Math.sqrt((block_sides_inner[0]+radius)*(block_sides_inner[0]+radius)+(block_sides_inner[1]+radius)*(block_sides_inner[1]+radius)+(block_sides_inner[2]+radius)*(block_sides_inner[2]+radius)); // Upper limit such that the pole plane is behind the box
      double[] xyz_pole1 = {ref_fr[0]+0.5*block_sides_outer[0]-shift*tilt[0], ref_fr[1]+0.5*block_sides_outer[1]-shift*tilt[1], ref_fr[2]+0.5*block_sides_outer[2]-shift*tilt[2]}; // Centre point on pole plane
      double[] vector_a = {Math.sin((theta[0]-90)*Math.PI/180)*Math.cos(theta[1]*Math.PI/180), Math.sin((theta[0]-90)*Math.PI/180)*Math.sin(theta[1]*Math.PI/180), Math.sqrt(1-Math.sin((theta[0]-90)*Math.PI/180)*Math.sin((theta[0]-90)*Math.PI/180))}; // Used for iteration in adda
      double[] vector_b = {tilt[1]*vector_a[2]-tilt[2]*vector_a[1], tilt[2]*vector_a[0]-tilt[0]*vector_a[2], tilt[0]*vector_a[1]-tilt[1]*vector_a[0]}; // Cross product of tilt and a
      int counter = 0;
      //fix - can improve performance by using the vectors to calc max adda and addb rather than using the limit
      double longest_length = 0.5*Math.sqrt(block_sides_outer[0]*block_sides_outer[0]+block_sides_outer[1]*block_sides_outer[1]+block_sides_outer[2]*block_sides_outer[2]);
      double adda = (longest_length+radius)/(3*radius);
      double addb = (longest_length+radius)/(3*radius);
      adda = Math.floor(adda);
      addb = Math.floor(addb);
      int addai = (int) adda;
      int addbi = (int) addb;
      int la = -addai;
      int ua = addai;
      int lb = -addbi;
      int ub = addbi;
      String[] allpoles = new String[(ua-la+1)*(ub-lb+1)];
      for (int a = la; a <= ua; a++) {
        for (int b = lb; b <= ub; b++) {
          //fix b
          double[] xyz_pole = {xyz_pole1[0]+(a*vector_a[0]+b*vector_b[0])*3*radius, xyz_pole1[1]+(a*vector_a[1]+b*vector_b[1])*3*radius, xyz_pole1[2]+(a*vector_a[2]+b*vector_b[2])*3*radius}; // Point of pole created using the original pole
          String poleTag = String.format("%s_pole_ax%d_ay%d", tag, a, b);
          model.geom("geom1").feature().create(poleTag, "Cylinder");
          model.geom("geom1").feature(poleTag).set("pos", xyz_pole); // Position of the cylinder
          model.geom("geom1").feature(poleTag).set("r", radius); // Cylinder radius
          model.geom("geom1").feature(poleTag).set("h", 2*shift); // Cylinder height
          model.geom("geom1").feature(poleTag).set("axis", tilt); // Axis
          allpoles[counter] = poleTag;
          counter++;
        }
      }
      
      
      
      // Make them a union
      String[] polesfilled = new String[counter];
      for (int n = 0; n < counter; n++) {
        polesfilled[n] = allpoles[n];
      }
      
      
      
      /// 3 Intersection
      String poleunionTag = tag+"_poles_union";
      model.geom("geom1").feature().create(poleunionTag, "Union");
      model.geom("geom1").feature(poleunionTag).selection("input").set(polesfilled);
      model.geom("geom1").feature(poleunionTag).set("intbnd", false); // Avoids unnecessary boundaries
      String intersectionTag = tag+"_intersection";
      model.geom("geom1").feature().create(intersectionTag, "Intersection");
      model.geom("geom1").feature(intersectionTag).selection("input").set(tempinnerTag, poleunionTag);
      
      
      
      /// 4 Build outside block
      String outerTag = tag+"_outer";
      model.geom("geom1").feature().create(outerTag, "Block");
      model.geom("geom1").feature(outerTag).set("base", "corner");
      model.geom("geom1").feature(outerTag).set("size", block_sides_outer);
      model.geom("geom1").feature(outerTag).set("pos", ref_fr);
      
      
      
      /// 5 Build inside block
      String innerTag = tag+"_inner";
      model.geom("geom1").feature().create(innerTag, "Block");
      model.geom("geom1").feature(innerTag).set("base", "corner");
      model.geom("geom1").feature(innerTag).set("size", block_sides_inner);
      model.geom("geom1").feature(innerTag).set("pos", block_inner);
      
      
      
      /// 6 Difference
      String differenceTag = tag+"_difference";
      model.geom("geom1").feature().create(differenceTag, "Difference");
      model.geom("geom1").feature(differenceTag).selection("input").set(outerTag);
      model.geom("geom1").feature(differenceTag).selection("input2").set(innerTag);
    }
  }
}



/// Generate Geometry
model.geom("geom1").run();



/// Add material
model.material().create("mat1", "Common", "mod1");
model.material("mat1").label("Matrix");
with(model.material("mat1").propertyGroup("def"));
  set("thermalconductivity", new String[]{"238[W/(m*K)]"}); set("density", new String[]{"2700[kg/m^3]"}); set("heatcapacity", new String[]{"900[J/(kg*K)]"});
endwith();
model.material("mat1").selection().all();



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
    model.selection(sel_coldsurface).set("xmax", xli);
    model.selection(sel_coldsurface).set("ymin", yli);
    model.selection(sel_coldsurface).set("ymax", yli);
    model.selection(sel_coldsurface).set("zmin", zli-0.01*wall_thickness);
    model.selection(sel_coldsurface).set("zmax", zli);
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