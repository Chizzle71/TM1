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



n = Integer.parseInt(model.param().get("n"));



/// Define parameters
String stdn = ""+n; // Number of the new study
//lattice sizing
int nucx = 12; //number of unit cells in x
int nucy = 12; //number of unit cells in y
int nucz = 12; //number of unit cells in z
double[] block_sides_outer = {0.08, 0.08, 0.08}; //size of unit cells (x,y,z)
double avgside = (block_sides_outer[0]+block_sides_outer[1]+block_sides_outer[2])/3;
//unit cell settings
double radius = 0.0012; //radius of poles
double rns = radius/avgside;
double wall_thickness = radius; //thickenss of walls
//inner block stuff
double[] theta = {0, 0};


////variables and building empty arrays
int nslices = nucz;
int n13 = 0;
String[][] inputthetas = new String[nslices][2];
double[][][] finput1 = new double[nucz][nucy][nucx];
double[][][] finput2 = new double[nucz][nucy][nucx];

////make array of the string names of the input.txt files
for (int n12 = 0; n12 < nucz; n12++) {
  n13 = n12+1+(n-1)*nucz;
  inputthetas[n12][0] = "input1thetas"+n13+".txt";
  inputthetas[n12][1] = "input2thetas"+n13+".txt";
}

////make lots of arrays in comsol of the actual data inside the input.txt files - this is a 12 by 12 - there are 2 * 12 of these
double[][] slice1input1 = readMatrixFromFile(inputthetas[0][0]);
double[][] slice1input2 = readMatrixFromFile(inputthetas[0][1]);

double[][] slice2input1 = readMatrixFromFile(inputthetas[1][0]);
double[][] slice2input2 = readMatrixFromFile(inputthetas[1][1]);

double[][] slice3input1 = readMatrixFromFile(inputthetas[2][0]);
double[][] slice3input2 = readMatrixFromFile(inputthetas[2][1]);

double[][] slice4input1 = readMatrixFromFile(inputthetas[3][0]);
double[][] slice4input2 = readMatrixFromFile(inputthetas[3][1]);

double[][] slice5input1 = readMatrixFromFile(inputthetas[4][0]);
double[][] slice5input2 = readMatrixFromFile(inputthetas[4][1]);

double[][] slice6input1 = readMatrixFromFile(inputthetas[5][0]);
double[][] slice6input2 = readMatrixFromFile(inputthetas[5][1]);

double[][] slice7input1 = readMatrixFromFile(inputthetas[6][0]);
double[][] slice7input2 = readMatrixFromFile(inputthetas[6][1]);

double[][] slice8input1 = readMatrixFromFile(inputthetas[7][0]);
double[][] slice8input2 = readMatrixFromFile(inputthetas[7][1]);

double[][] slice9input1 = readMatrixFromFile(inputthetas[8][0]);
double[][] slice9input2 = readMatrixFromFile(inputthetas[8][1]);

double[][] slice10input1 = readMatrixFromFile(inputthetas[9][0]);
double[][] slice10input2 = readMatrixFromFile(inputthetas[9][1]);

double[][] slice11input1 = readMatrixFromFile(inputthetas[10][0]);
double[][] slice11input2 = readMatrixFromFile(inputthetas[10][1]);

double[][] slice12input1 = readMatrixFromFile(inputthetas[11][0]);
double[][] slice12input2 = readMatrixFromFile(inputthetas[11][1]);



/////takes variables from 2 dimension arrays and puts them into a 3 dimension array
double tempav1 = 0;
double tempav2 = 0;

for (int j = 0; j < nucy; j++) {
  for (int k = 0; k < nucx; k++) {
    tempav1 = slice1input1[j][k];
    tempav2 = slice1input2[j][k];
    finput1[0][j][k] = tempav1;
    finput2[0][j][k] = tempav2;
    
    tempav1 = slice2input1[j][k];
    tempav2 = slice2input2[j][k];
    finput1[1][j][k] = tempav1;
    finput2[1][j][k] = tempav2;
    
    tempav1 = slice3input1[j][k];
    tempav2 = slice3input2[j][k];
    finput1[2][j][k] = tempav1;
    finput2[2][j][k] = tempav2;
    
    tempav1 = slice4input1[j][k];
    tempav2 = slice4input2[j][k];
    finput1[3][j][k] = tempav1;
    finput2[3][j][k] = tempav2;
    
    tempav1 = slice5input1[j][k];
    tempav2 = slice5input2[j][k];
    finput1[4][j][k] = tempav1;
    finput2[4][j][k] = tempav2;
    
    tempav1 = slice6input1[j][k];
    tempav2 = slice6input2[j][k];
    finput1[5][j][k] = tempav1;
    finput2[5][j][k] = tempav2;
    
    tempav1 = slice7input1[j][k];
    tempav2 = slice7input2[j][k];
    finput1[6][j][k] = tempav1;
    finput2[6][j][k] = tempav2;
    
    tempav1 = slice8input1[j][k];
    tempav2 = slice8input2[j][k];
    finput1[7][j][k] = tempav1;
    finput2[7][j][k] = tempav2;
    
    tempav1 = slice9input1[j][k];
    tempav2 = slice9input2[j][k];
    finput1[8][j][k] = tempav1;
    finput2[8][j][k] = tempav2;
    
    tempav1 = slice10input1[j][k];
    tempav2 = slice10input2[j][k];
    finput1[9][j][k] = tempav1;
    finput2[9][j][k] = tempav2;
    
    tempav1 = slice11input1[j][k];
    tempav2 = slice11input2[j][k];
    finput1[10][j][k] = tempav1;
    finput2[10][j][k] = tempav2;
    
    tempav1 = slice12input1[j][k];
    tempav2 = slice12input2[j][k];
    finput1[11][j][k] = tempav1;
    finput2[11][j][k] = tempav2;
  }
}



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



/// Loop through block
for (int i = 0; i < nucz; i++) {
  for (int j = 0; j < nucy; j++) {
    for (int k = 0; k < nucx; k++) {
      ////the final thing - takes the values from the two final input arrays for the exact unit cell it is in inside the loop
      theta[0] = finput1[i][j][k];
      theta[1] = finput2[i][j][k];
      
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



/// Results
double plane_step = block_sides_outer[2]/4;
double nplanes = (nucz*block_sides_outer[2])/(plane_step)+1;
double z_initial = 0;
String tagd = "datad"+n;
String tage = "datae"+n;
for (int i = 0; i < nplanes; i++) {
  String tagdi = tagd+0.001*i;
  String tagei = tage+0.001*i;
  double z_value = z_initial+i*plane_step;
  model.result().dataset().create(tagdi, "CutPlane"); //creates the dataset plane
  model.result().dataset(tagdi).set("data", "dset"+n); //sets data to take from datasetn (solution)
  model.result().dataset(tagdi).set("quickplane", "xy"); //set dataset to xy plane
  model.result().dataset(tagdi).set("quickz", z_value); //set z value
  model.result().export().create(tagei, "Data"); //create export
  model.result().export(tagei).set("data", tagdi); //set data to tagd which is the plane just created
  model.result().export(tagei).set("filename", "output"+n+"level"+i+".txt"); //set filename
  model.result().export(tagei).run(); //run
}