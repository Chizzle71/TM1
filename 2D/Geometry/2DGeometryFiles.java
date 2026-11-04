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


String file = model.param().get("inputfile");
int n = (int) model.param().evaluate("n");


// Define parameters
String stdn = ""+n; // Number of the new study
int nucx = 2;
int nucy = 2;
double sidex = 0.01; // Rectangle length in x
double sidey = 0.01; // Rectangle length in y
double sidem = sidey;
if (sidex > sidey) {
  sidem = sidex;
}
double h = 0.0004; // Cut out width
double h_y = 0;
double h_x = 0;
double C_i = 0;
double D_i = 0;
double cbp = sidem/h;
int arraysize = (int) cbp+1;
double[][] As = new double[arraysize][2];
double[][] Bs = new double[arraysize][2];
double[][] Cs = new double[arraysize][2];
double[][] Ds = new double[arraysize][2];
double[] xs = {0, 0, 0, 0, 0, 0};
double[] ys = {0, 0, 0, 0, 0, 0};
String[] innerTags = new String[0];
double linecount = nucx*(nucy+1)+nucy*(nucx+1);


double[][] display_grid = readMatrixFromFile(file);


// Loop through grid and create squares
for (int j = 0; j < nucy; j++) {
  for (int i = 0; i < nucx; i++) {
    double theta = display_grid[j][i];
    String tag = String.format("sq_%d_%d", i, j);
    
    if (theta == 360) {
      continue;
    }
    
    // Unit cell
    // Variables
    double x = i*sidex; // Spacing
    double y = j*sidey; // Spacing
    int An = 0; // Polypoint count
    int Bn = 0; // Polypoint count
    int Cn = 0; // Polypoint count
    int Dn = 0; // Polypoint count
    int m = 1; // Polarity of theta
    int ignore = 0;
    
    // Check for negative theta
    if (theta < 0) {
      theta = -theta;
      m = -m;
    }
    
    // Define second level parameters
    if (theta == 90) {
      h_y = sidey;
      h_x = h;
      ignore = 1;
    }
    else if (theta == 0) {
      h_y = h;
      h_x = sidex;
      ignore = 1;
    }
    else {
      h_y = h/Math.cos(theta*Math.PI/180);
      h_x = h/Math.sin(theta*Math.PI/180);
    }
    
    // Make outer square
    String outerTag = tag+"_outer";
    model.geom("geom1").feature().create(outerTag, "Rectangle");
    model.geom("geom1").feature(outerTag).set("base", "corner");
    model.geom("geom1").feature(outerTag).set("size", new double[]{sidex, sidey});
    model.geom("geom1").feature(outerTag).set("pos", new double[]{x, y});
    
    // Define polygon points
    
    // Left side
    for (double A = sidey-h; A >= h; A = A-h_y) {
      // Add this point to the array
      As[An][0] = h;
      As[An][1] = A;
      // Increase polypoint count
      An++;
      // Variables for later
      C_i = h+h_x*(1-(A-h)/(h_y));
    }
    
    // Top side
    for (double B = h; B <= sidex-h; B = B+h_x) {
      // Add this point to the array
      Bs[Bn][0] = B;
      Bs[Bn][1] = sidey-h;
      // Increase polypoint count
      Bn++;
      // Variables for later
      if (Bn == 1) {
        D_i = sidey-2*h+(sidex-2*h)*Math.tan(theta*Math.PI/180);
      }
      else if (Bn > 1) {
        D_i = sidey-h-h_y*(1-(sidex-h-B)/(h_x));
      }
    }
    
    // Bottom side
    if (theta == 90) {C_i = 2*h; }
    for (double C = C_i; C <= sidex-h; C = C+h_x) {
      // Add this point to the array
      Cs[Cn][0] = C;
      Cs[Cn][1] = h;
      // Increase polypoint count
      Cn++;
    }
    
    // Right side
    for (double D = D_i; D >= h; D = D-h_y) {
      // Add this point to the array
      Ds[Dn][0] = sidex-h;
      Ds[Dn][1] = D;
      // Increase polypoint count
      Dn++;
    }
    
    // Make polygon geometry array
    double npolys = ((sidex-2*h)/h_x+(sidey-2*h)/h_y)/2-ignore;
    linecount = linecount+4*npolys-1;
    for (int c2 = 0; c2 <= npolys; c2++) {
      
      // For points along A and C
      if (2*c2+1-An < 0) {
        // 2 points on A
        xs[0] = h;
        xs[1] = h;
        xs[2] = h;
        ys[0] = As[2*c2][1];
        ys[1] = As[2*c2][1];
        ys[2] = As[2*c2+1][1];
      }
      else if (2*c2+1-An == 0) {
        // 1 point on A, 1 point on C + bottom left corner
        xs[0] = h;
        xs[1] = h;
        xs[2] = Cs[0][0];
        ys[0] = As[2*c2][1];
        ys[1] = h;
        ys[2] = h;
        linecount++;
      }
      else {
        // Only points on C
        if (c2*2-An-Cn+1 < 0) {
          // No corner point
          xs[0] = Cs[c2*2-An][0];
          xs[1] = Cs[c2*2-An+1][0];
          xs[2] = Cs[c2*2-An+1][0];
          ys[0] = h;
          ys[1] = h;
          ys[2] = h;
        }
        else {
          // With corner point
          xs[0] = Cs[c2*2-An][0];
          xs[1] = sidex-h;
          xs[2] = sidex-h;
          ys[0] = h;
          ys[1] = h;
          ys[2] = h;
          linecount--;
        }
      }
      
      // For points along B and D
      if (2*c2+1-Bn < 0) {
        // 2 points on B
        xs[5] = Bs[c2*2][0];
        xs[4] = Bs[c2*2][0];
        xs[3] = Bs[c2*2+1][0];
        ys[4] = sidey-h;
        ys[5] = sidey-h;
        ys[3] = sidey-h;
      }
      else if (2*c2+1-Bn == 0) {
        // 1 point on B, 1 point on D + top right corner
        xs[5] = Bs[c2*2][0];
        xs[4] = sidex-h;
        xs[3] = sidex-h;
        ys[5] = sidey-h;
        ys[4] = sidey-h;
        ys[3] = Ds[0][1];
        linecount++;
      }
      else {
        // Only points on D
        if (c2*2-Bn-Dn < 0) {
          // No corner point
          xs[5] = sidex-h;
          xs[4] = sidex-h;
          xs[3] = sidex-h;
          ys[5] = Ds[c2*2-Bn][1];
          ys[4] = Ds[c2*2-Bn+1][1];
          ys[3] = Ds[c2*2-Bn+1][1];
        }
        else {
          // With corner point
          xs[5] = sidex-h;
          xs[4] = sidex-h;
          xs[3] = sidex-h;
          ys[5] = Ds[c2*2-Bn][1];
          ys[4] = h;
          ys[3] = h;
        }
      }
      
      // Add x and y to make change frame from relative to absolute frame. And apply symmetry if negative
      if (m == 1) {
        for (int m1 = 0; m1 < 6; m1++) {
          xs[m1] = x+xs[m1];
          ys[m1] = y+ys[m1];
        }
      }
      else if (m == -1) {
        for (int m1 = 0; m1 < 6; m1++) {
          xs[m1] = x+sidex-xs[m1];
          ys[m1] = y+ys[m1];
        }
      }
      
      // Make polygon as solid cutout
      String innerTag = "cutout_"+i+"_"+j+"_"+c2;
      String[] newArr = new String[innerTags.length+1];
      for (int k = 0; k < innerTags.length; k++) {
        newArr[k] = innerTags[k];
      }
      newArr[innerTags.length] = innerTag;
      innerTags = newArr;
      model.geom("geom1").feature().create(innerTag, "Polygon"); //geometry type 1
      model.geom("geom1").feature(innerTag).set("x", xs); //set x values
      model.geom("geom1").feature(innerTag).set("y", ys); //set y values
      model.geom("geom1").feature(innerTag).set("type", "solid"); //set geometry type 2
    }
    
    // Difference = square - polygons
    String diffTag = tag+"_diff";
    model.geom("geom1").feature().create(diffTag, "Difference");
    model.geom("geom1").feature(diffTag).selection("input").set(new String[]{outerTag});
    model.geom("geom1").feature(diffTag).selection("input2").set(innerTags);
    
  }
}

// Add geometry
model.geom("geom1").run();


// Add material
model.material().create("mat1", "Common", "mod1");
model.material("mat1").label("Matrix");
with(model.material("mat1").propertyGroup("def"));
  set("thermalconductivity", new String[]{"238[W/(m*K)]"}); set("density", new String[]{"2700[kg/m^3]"}); set("heatcapacity", new String[]{"900[J/(kg*K)]"});
endwith();
model.material("mat1").selection().all();


// Calculations for physics
for (int i = 0; i < nucy; i++) {
  String sel_hotline = "hline"+i;
  double yli = i*sidey+h;
  model.selection().create(sel_hotline, "Box");
  model.selection(sel_hotline).geom("geom1", 1);
  model.selection(sel_hotline).set("xmin", -0.01*h);
  model.selection(sel_hotline).set("xmax", 0.01*h);
  model.selection(sel_hotline).set("ymin", yli);
  model.selection(sel_hotline).set("ymax", yli+h);
}
for (int i = 0; i < nucy; i++) {
  String sel_coldline = "cline"+i;
  double xli = nucx*sidex;
  double yli = i*sidey+h;
  model.selection().create(sel_coldline, "Box");
  model.selection(sel_coldline).geom("geom1", 1);
  model.selection(sel_coldline).set("xmin", xli-0.01*h);
  model.selection(sel_coldline).set("xmax", xli+h);
  model.selection(sel_coldline).set("ymin", yli);
  model.selection(sel_coldline).set("ymax", yli+h);
}


// Add physics
model.physics().create("ht", "HeatTransfer", "geom1");
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
  String temp = "ctemp"+j;
  model.physics("ht").create(temp, "TemperatureBoundary", 1);
  model.physics("ht").feature(temp).selection().named(sel_coldline);
  model.physics("ht").feature(temp).set("T0", "300[K]");
}


// Add mesh
model.mesh().create("mesh1", "geom1");
model.mesh("mesh1").automatic(true); // Keep on Physics-controlled mesh
model.mesh("mesh1").autoMeshSize(3); // [1,9] [fine,coarse]
model.mesh("mesh1").run();


// Add study
model.study().create(stdn);
model.study(stdn).create("stat", "Stationary");
model.study(stdn).run();


// Results
model.result().export().create("data"+n, "Data");
model.result().export("data"+n).set("data", "dset"+n);
model.result().export("data"+n).set("filename", "output"+n+".txt");
model.result().export("data"+n).run();