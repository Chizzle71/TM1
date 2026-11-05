// Remove old geometry if it exists
String[] geomTags = model.geom().tags();
for (String tag : geomTags) {
  if (tag.equals("geom1")) {
    model.geom().remove(tag);
    break;
  }
}


// Create new 2D geometry
model.geom().create("geom1", 2);
model.geom("geom1").lengthUnit("m");


// Define parameters
String stdn = "4"; // Number of the new study + 1
int nucx = 30;
int nucy = 30;
double side = 0.04; // Square length
double h = 0.002; // Cut out width
double E = 0.218; // Efficiency guess


//hand
double[][] display_grid = {
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
  {45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45}
};


// Loop through grid and create squares
for (int j = 0; j < nucy; j++) {
  for (int i = 0; i < nucx; i++) {
    double theta = display_grid[j][i];
    String tag = String.format("sq_%d_%d", i, j);
   
    // Unit cell
    // Variables
    double x = i*side; // Spacing
    double y = j*side; // Spacing
   
    if (theta == 360) {
      continue;
    }
   
    // Make outer square
    String outerTag = tag+"_outer";
    model.geom("geom1").feature().create(outerTag, "Square");
    model.geom("geom1").feature(outerTag).set("base", "corner");
    model.geom("geom1").feature(outerTag).set("size", side);
    model.geom("geom1").feature(outerTag).set("pos", new double[]{x, y});
   
    // Calculate effective conductivity matrix
    double thetar = theta*Math.PI/180;
    double kxx = 238*(Math.cos(thetar)*Math.cos(thetar)+E*Math.sin(thetar)*Math.sin(thetar));
    double kyx = 238*((1-E)*Math.cos(thetar)*Math.sin(thetar));
    double kxy = kyx; // Symmetry condition on
    double kyy = 238*(E*Math.cos(thetar)*Math.cos(thetar)+Math.sin(thetar)*Math.sin(thetar));
   
    // Apply conductivity matrix
    //    model.physics().create("ht"+tag, "HeatTransfer", "geom1");
    //    String solid = "solid"+i+j*0.001;
    //    model.physics("ht").create("sdd", "Domain");
    //    model.physics("ht").feature(solid).set("k_mat", "userdef");
   
    String sm = "sm"+i+0000+j;
    model.selection().create(sm, "Box");
    model.selection(sm).geom("geom1", 2);
    model.selection(sm).set("xmin", x+h);
    model.selection(sm).set("xmax", x+h*2);
    model.selection(sm).set("ymin", y+h);
    model.selection(sm).set("ymax", y+h*2);
   
    String mat = "mat"+i+0000+j;
    model.material().create(mat, "Common", "mod1");
    model.material(mat).label(mat);
    with(model.material(mat).propertyGroup("def"));
      set("thermalconductivity", new String[]{
        kxx+"[W/(m*K)]", kyx+"[W/(m*K)]", "0",
        kxy+"[W/(m*K)]", kyy+"[W/(m*K)]", "0",
        "0", "0", "0"
      });
      set("density", new String[]{"2700[kg/m^3]"});
      set("heatcapacity", new String[]{"900[J/(kg*K)]"});
    endwith();
    model.material(mat).selection().named(sm);
   
  }
}


// Add geometry
model.geom("geom1").run();


// Physics
model.physics().create("ht", "HeatTransfer", "geom1");
//model.physics("ht").feature("solid1").set("k_mat", "userdef");
//model.physics("ht").feature("solid1").set("k", new double[]{
//  1, 0, 0.0, // kxx kxy kxz
//  0, 0.218, 0.0, // kyx kyy kyz
//  0.0, 0.0, 0.0 // kzx kzy kzz
//});
//calculations
for (int i = 0; i < nucy; i++) {
  String sel_hotline = "hline"+i;
  double yli = i*side+h;
  model.selection().create(sel_hotline, "Box");
  model.selection(sel_hotline).geom("geom1", 1);
  model.selection(sel_hotline).set("xmin", -0.01*h);
  model.selection(sel_hotline).set("xmax", 0.01*h);
  model.selection(sel_hotline).set("ymin", yli);
  model.selection(sel_hotline).set("ymax", yli+h);
}
for (int i = 0; i < nucy; i++) {
  String sel_coldline = "cline"+i;
  double xli = nucx*side;
  double yli = i*side+h;
  model.selection().create(sel_coldline, "Box");
  model.selection(sel_coldline).geom("geom1", 1);
  model.selection(sel_coldline).set("xmin", xli-0.01*h);
  model.selection(sel_coldline).set("xmax", xli+h);
  model.selection(sel_coldline).set("ymin", yli);
  model.selection(sel_coldline).set("ymax", yli+h);
}
//hot
for (int i = 0; i < nucy; i++) {
  int j = i+2;
  String sel_hotline = "hline"+i;
  String temp = "htemp"+j;
  model.physics("ht").create(temp, "TemperatureBoundary", 1);
  model.physics("ht").feature(temp).selection().named(sel_hotline);
  model.physics("ht").feature(temp).set("T0", "400[K]");
}
//cold
for (int i = 0; i < nucy; i++) {
  int j = i+2;
  String sel_coldline = "cline"+i;
  String temp = "temp"+j;
  model.physics("ht").create(temp, "TemperatureBoundary", 1);
  model.physics("ht").feature(temp).selection().named(sel_coldline);
  model.physics("ht").feature(temp).set("T0", "300[K]");
}


// Add mesh
model.mesh().create("mesh1", "geom1");
model.mesh("mesh1").automatic(true); // Keep on Physics-controlled mesh
model.mesh("mesh1").autoMeshSize(2); // [1,9] [fine,coarse]
model.mesh("mesh1").run();


// Add study
model.study().create(stdn);
model.study(stdn).create("stat", "Stationary");
model.study(stdn).run();