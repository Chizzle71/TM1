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
int nucx = 6;
int nucy = 8;
double side = 0.03; // Square length
double h = 0.0011; // Cut out width
double h_y = 0;
double h_x = 0;
double C_i = 0;
double D_i = 0;
double cbp = side/h;
int arraysize = (int) cbp;
double[][] As = new double[arraysize][2];
double[][] Bs = new double[arraysize][2];
double[][] Cs = new double[arraysize][2];
double[][] Ds = new double[arraysize][2];
double[] xs = {0, 0, 0, 0, 0, 0};
double[] ys = {0, 0, 0, 0, 0, 0};
String[] innerTags = new String[0];

//double[][] display_grid = {
//  {-45, -45, -45, -45, -45, -45},
//  {-45, -45, -45, -45, -45, -45},
//  {-45, -45, -45, -45, -45, -45},
//  {-45, -45, -45, -45, -45, -45},
//  {45, 45, 45, 45, 45, 45},
//  {45, 45, 45, 45, 45, 45},
//  {45, 45, 45, 45, 45, 45},
//  {45, 45, 45, 45, 45, 45}
//};

double[][] grid = new double[nucxy][nucxy];
for (int i = 0; i < nucxy; i++) {
  for (int j = 0; j < nucxy; j++) {
    grid[i][j] = display_grid[nucxy-1-j][i]; // Matrix transformation
  }
}

// Loop through grid and create squares
for (int j = 0; j < nucy; j++) {
  for (int i = 0; i < nucx; i++) {
    double theta = display_grid[j][i];
    String tag = String.format("sq_%d_%d", i, j);
    
    // Unit cell
    // Variables
    double x = i*side; // Spacing
    double y = j*side; // Spacing
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
      h_y = side;
      h_x = h;
      ignore = 1;
    }
    else if (theta == 0) {
      h_y = h;
      h_x = side;
      ignore = 1;
    }
    else {
      h_y = h/Math.cos(theta*Math.PI/180);
      h_x = h/Math.sin(theta*Math.PI/180);
    }
    
    // Make outer square
    String outerTag = tag+"_outer";
    model.geom("geom1").feature().create(outerTag, "Square");
    model.geom("geom1").feature(outerTag).set("base", "corner");
    model.geom("geom1").feature(outerTag).set("size", side);
    model.geom("geom1").feature(outerTag).set("pos", new double[]{x, y});
    
    // Define polygon points
    
    // Left side
    for (double A = side-h; A >= h; A = A-h_y) {
      // Add this point to the array
      As[An][0] = h;
      As[An][1] = A;
      // Increase polypoint count
      An++;
      // Variables for later
      C_i = h+h_x*(1-(A-h)/(h_y));
    }
    
    // Top side
    for (double B = h; B <= side-h; B = B+h_x) {
      // Add this point to the array
      Bs[Bn][0] = B;
      Bs[Bn][1] = side-h;
      // Increase polypoint count
      Bn++;
      // Variables for later
      if (Bn == 1) {
        D_i = side-2*h+(side-2*h)*Math.tan(theta*Math.PI/180);
      }
      else if (Bn > 1) {
        D_i = side-h-h_y*(1-(side-h-B)/(h_x));
      }
    }
    
    // Bottom side
    if (theta == 90) {C_i = 2*h; }
    for (double C = C_i; C <= side-h; C = C+h_x) {
      // Add this point to the array
      Cs[Cn][0] = C;
      Cs[Cn][1] = h;
      // Increase polypoint count
      Cn++;
    }
    
    // Right side
    for (double D = D_i; D >= h; D = D-h_y) {
      // Add this point to the array
      Ds[Dn][0] = side-h;
      Ds[Dn][1] = D;
      // Increase polypoint count
      Dn++;
    }
    
    // Make polygon geometry array
    double npolys = ((side-2*h)/h_x+(side-2*h)/h_y)/2-ignore;
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
          xs[1] = side-h;
          xs[2] = side-h;
          ys[0] = h;
          ys[1] = h;
          ys[2] = h;
        }
      }
      
      // For points along B and D
      if (2*c2+1-Bn < 0) {
        // 2 points on B
        xs[5] = Bs[c2*2][0];
        xs[4] = Bs[c2*2][0];
        xs[3] = Bs[c2*2+1][0];
        ys[4] = side-h;
        ys[5] = side-h;
        ys[3] = side-h;
      }
      else if (2*c2+1-Bn == 0) {
        // 1 point on B, 1 point on D + top right corner
        xs[5] = Bs[c2*2][0];
        xs[4] = side-h;
        xs[3] = side-h;
        ys[5] = side-h;
        ys[4] = side-h;
        ys[3] = Ds[0][1];
      }
      else {
        // Only points on D
        if (c2*2-Bn-Dn < 0) {
          // No corner point
          xs[5] = side-h;
          xs[4] = side-h;
          xs[3] = side-h;
          ys[5] = Ds[c2*2-Bn][1];
          ys[4] = Ds[c2*2-Bn+1][1];
          ys[3] = Ds[c2*2-Bn+1][1];
        }
        else {
          // With corner point
          xs[5] = side-h;
          xs[4] = side-h;
          xs[3] = side-h;
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
          xs[m1] = x+side-xs[m1];
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

// Finalise geometry
model.geom("geom1").run();

// Add material
// Get the material node (create it if it doesn't exist)
model.material().create("mat1", "Common");
model.material("mat1").label("Copper");

// Set the material from the COMSOL Material Library
model.material("mat1").set("family", "copper");

// Assign material to all domains
model.material("mat1").selection().all();


// Add mesh
model.mesh().create("mesh1", "geom1");
model.mesh("mesh1").automatic(true); // Turn off Physics-controlled mesh with(model.mesh("mesh1").feature("size"));
model.mesh("mesh1").autoMeshSize(7); // [1,9] [fine,coarse]
model.mesh("mesh1").run();