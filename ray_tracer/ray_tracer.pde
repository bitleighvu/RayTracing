// This is the starter code for the CS 3451 Ray Tracing Project.
//
// The most important part of this code is the interpreter, which will help
// you parse the scene description (.cli) files.

// BitLeigh Vu

boolean debug_flag;  // help to print information for just one pixel

ArrayList<Light> lights = new ArrayList<Light>();
ArrayList<Disk> disks = new ArrayList<Disk>();
ArrayList<Ellipsoid> ellipsoids = new ArrayList<Ellipsoid>();
Surface[] sur = new Surface[5];

PVector background;
PVector eye;
PVector u;
PVector v;
PVector w;
float d;

void setup() {
  size(500, 500);  
  noStroke();
  colorMode(RGB);
  background(0, 0, 0);
}

void reset_scene() {
  //reset the global scene variables here
  sur = new Surface[5];
  lights = new ArrayList<Light>();
  disks = new ArrayList<Disk>();
  ellipsoids = new ArrayList<Ellipsoid>();
  d = 0;
  background = new PVector(0, 0, 0);
  eye = new PVector(0, 0, 0);
  u = new PVector(0, 0, 0);
  v = new PVector(0, 0, 0);
  w = new PVector(0, 0, 0);
}

void keyPressed() {
  reset_scene();
  switch(key) {
  case '1':  
    interpreter("01.cli"); 
    break;
  case '2':  
    interpreter("02.cli"); 
    break;
  case '3':  
    interpreter("03.cli"); 
    break;
  case '4':  
    interpreter("04.cli"); 
    break;
  case '5':  
    interpreter("05.cli"); 
    break;
  case '6':  
    interpreter("06.cli"); 
    break;
  case '7':  
    interpreter("07.cli"); 
    break;
  case '8':  
    interpreter("08.cli"); 
    break;
  case '9':  
    interpreter("09.cli"); 
    break;
  case '0':  
    interpreter("10.cli"); 
    break;
  case '-':  
    interpreter("11.cli"); 
    break;
  case 'q':  
    exit(); 
    break;
  }
}

// this routine helps parse the text in the scene description files
void interpreter(String filename) {

  println("Parsing '" + filename + "'");
  String str[] = loadStrings(filename);
  if (str == null) println("Error! Failed to read the cli file.");

  for (int i = 0; i < str.length; i++) {

    String[] token = splitTokens(str[i], " ");  // Get a line and parse the tokens

    if (token.length == 0) continue; // Skip blank lines

    if (token[0].equals("fov")) {
      float fov = float(token[1]);

      // call routine to save the fov
      d = 1 / (tan((radians(fov))/2));
    } else if (token[0].equals("background")) {
      float r = float(token[1]);
      float g = float(token[2]);
      float b = float(token[3]);

      // call routine to save the background color
      background = new PVector(r, g, b);
    } else if (token[0].equals("eye")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);

      // call routine to save the eye position
      eye = new PVector(x, y, z);
    } else if (token[0].equals("uvw")) {
      float ux = float(token[1]);
      float uy = float(token[2]);
      float uz = float(token[3]);
      float vx = float(token[4]);
      float vy = float(token[5]);
      float vz = float(token[6]);
      float wx = float(token[7]);
      float wy = float(token[8]);
      float wz = float(token[9]);

      // call routine to save the camera's values for u,v,w
      u = new PVector(ux, uy, uz);
      v = new PVector(vx, vy, vz);
      w = new PVector(wx, wy, wz);
    } else if (token[0].equals("light")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);
      float r = float(token[4]);
      float g = float(token[5]);
      float b = float(token[6]);

      // call routine to save lighting information
      Light light = new Light(x, y, z, r, g, b);
      lights.add(light);
    } else if (token[0].equals("surface")) {
      float Cdr = float(token[1]);
      float Cdg = float(token[2]);
      float Cdb = float(token[3]);
      float Car = float(token[4]);
      float Cag = float(token[5]);
      float Cab = float(token[6]);
      float Csr = float(token[7]);
      float Csg = float(token[8]);
      float Csb = float(token[9]);
      float P = float(token[10]);
      float K = float(token[11]);


      // call routine to save the surface material properties
      Surface surface = new Surface(Cdr, Cdg, Cdb, Car, Cag, Cab, Csr, Csg, Csb, P, K);
      sur[1] = surface;
    } else if (token[0].equals("ellipsoid")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);
      float rx = float(token[4]);
      float ry = float(token[5]);
      float rz = float(token[6]);

      // call routine to save ellipsoid here
      Ellipsoid ellipsoid = new Ellipsoid(x, y, z, rx, ry, rz, sur[1]);
      ellipsoids.add(ellipsoid);
    } else if (token[0].equals("disk")) {
      float x = float(token[1]);
      float y = float(token[2]);
      float z = float(token[3]);
      float nx = float(token[4]);
      float ny = float(token[5]);
      float nz = float(token[6]);
      float radius = float(token[7]);

      // call routine to save disk here
      Disk disk = new Disk(x, y, z, nx, ny, nz, radius, sur[1]);
      disks.add(disk);
    } else if (token[0].equals("write")) {
      draw_scene();   // here is where you actually perform the ray tracing
      println("Saving image to '" + token[1] + "'");
      save(token[1]); // this saves your ray traced scene to a .png file
    } else if (token[0].equals("#")) {
      // comment symbol (ignore this line)
    } else {
      println ("cannot parse this line: " + str[i]);
    }
  }
}

// This is where you should place your code for creating the eye rays and tracing them.
void draw_scene() {  
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {

      // maybe turn on a debug flag for a particular pixel (so you can print ray information)
      if (x == 84 && y == 197)
        debug_flag = true;
      else
        debug_flag = false;

      // create and cast an eye ray here
      Ray eyeRay = createEyeRay(x, y);

      Hit hit = intersection(eyeRay);
      PVector finalColor = calcColor(hit, 1);
      fill(finalColor.x * 255.0, finalColor.y * 255.0, finalColor.z * 255.0);
      rect (x, y, 1, 1);
    }
  }
}

void draw() {
  // nothing should be put here for this project
}

// use this routine to find the coordinates of a particular pixel (for debugging)
void mousePressed() {
  println ("Mouse pressed at location: " + mouseX + " " + mouseY);
}

// calculates the color for hit point and recursive reflection
PVector calcColor(Hit hit, int depth) {
  if (hit == null) {
    return background;
  }

  PVector colors = hit.shape.ambient;

  for (Light light : lights) {
    // diffuse color
    PVector l = PVector.sub(light.pos, hit.intersect).normalize();
    float normalL = max(0, PVector.dot(hit.normal, l));
    PVector diff = PVector.mult(light.colors, normalL);
    
    // specular color
    PVector eyeRay = PVector.sub(hit.ray.origin, hit.intersect).normalize();
    float hDotN = pow(PVector.dot(PVector.add(eyeRay, l).normalize(), hit.normal), hit.shape.phongExp);
    PVector spec = PVector.mult(light.colors, hDotN); 
    
    // shadow ray 
    PVector origin = PVector.add(hit.intersect, PVector.mult(hit.normal, 0.0001));
    PVector direction = PVector.sub(light.pos, hit.intersect).normalize();
    
    if (intersection(new Ray(origin, direction)) == null || hit.tVal < 1) {
      colors = PVector.add(colors, new PVector(diff.x * hit.shape.diffuse.x, diff.y * hit.shape.diffuse.y, diff.z * hit.shape.diffuse.z));
      colors = PVector.add(colors, new PVector(spec.x * hit.shape.specular.x, spec.y * hit.shape.specular.y, spec.z * hit.shape.specular.z));  
    }
  }
  
  // reflection
  PVector reflection_color = new PVector(0, 0, 0);
  if (hit.shape.Krefl <= 1 && hit.shape.Krefl > 0 && depth < 10) {
    PVector rayOrigin = PVector.add(hit.intersect, PVector.mult(hit.normal, 0.0001));
    PVector eyeRay = PVector.sub(hit.ray.origin, hit.intersect).normalize();
    float NdotE = PVector.dot(hit.normal, eyeRay);
    PVector R = PVector.sub(PVector.mult(hit.normal, 2 * NdotE), eyeRay).normalize();
    
    reflection_color = calcColor(intersection(new Ray(rayOrigin, R)), depth + 1);
  }

  return PVector.add(colors, PVector.mult(reflection_color, hit.shape.Krefl));  
  
}

// Ray class
class Ray {
  PVector origin, direction;

  Ray(PVector eye, PVector direc) {
    origin = eye;
    direction = direc;
  }
}


Ray createEyeRay(int x, int y) {
  // ray direction = -dw + uu - vv
  float uPoint = ((2.0 * x)/ (float) width) - 1.0;
  float vPoint = ((2.0 * y)/ (float) height) - 1.0;

  PVector first = PVector.mult(w, -d);
  PVector second = PVector.mult(u, uPoint);
  PVector third = PVector.mult(v, vPoint); 

  PVector eyeDir = PVector.add(first, second);
  eyeDir = eyeDir.sub(third);

  Ray eyeRay = new Ray(eye, eyeDir); // position of Eye, eyedirection
  return eyeRay;
}

// hit info
class Hit {
  Shape shape;
  float tVal;
  PVector intersect, normal;
  Ray ray;

  Hit(Shape shapeP, float t, PVector inter, PVector norm, Ray eyeR) {
    shape = shapeP;
    tVal = t;
    intersect = inter;
    normal = norm;
    ray = eyeR;
  }
}

Hit intersection(Ray ray) {
  float t = Float.MAX_VALUE;
  float minRoot = Float.MAX_VALUE;
  Hit hit = null;

  // disk intersection
  for (Disk disk : disks) {
    // vector for projection onto plane
    PVector proj = PVector.sub(disk.pos, ray.origin); // disk center - eye origin
    float top = PVector.dot(proj, disk.normal); // (projection vector * plane normal)
    float bottom = PVector.dot(disk.normal, ray.direction); // (disk normal * eye direction)
    float curT = top / bottom;
    
    PVector intercept = new PVector(ray.origin.x + (ray.direction.x * curT), ray.origin.y + (ray.direction.y * curT), ray.origin.z + (ray.direction.z * curT));
    float dist = PVector.dist(intercept, disk.pos);

    if (curT < t && curT >= 0) {
      if (dist <= disk.radius) {
        t = curT;
        hit = new Hit(disk, t, intercept, disk.normal, ray); // Shape shapeP, float t, PVector inter, PVector norm, Ray eyeR
      }
    }
  }
 
  // ellipsoid intersection
  for (Ellipsoid ellipsoid : ellipsoids) {   
    // difference between where ray origin is and ellipsoid 
    PVector delta = PVector.sub(ray.origin, ellipsoid.pos);
    
     float a = (pow(ray.direction.x, 2) / pow(ellipsoid.radi.x, 2)) +
       (pow(ray.direction.y, 2) / pow(ellipsoid.radi.y, 2)) +
       (pow(ray.direction.z, 2) / pow(ellipsoid.radi.z, 2));
     
     float b = 2 * (
       ((delta.x * ray.direction.x) / pow(ellipsoid.radi.x, 2)) +
       ((delta.y * ray.direction.y) / pow(ellipsoid.radi.y, 2)) +
       ((delta.z * ray.direction.z) / pow(ellipsoid.radi.z, 2))
       );
       
    float c = (pow(delta.x, 2) / pow(ellipsoid.radi.x, 2)) +
      (pow(delta.y, 2) / pow(ellipsoid.radi.y, 2)) +
      (pow(delta.z, 2) / pow(ellipsoid.radi.z, 2)) 
      - 1;
  
    float determinant = pow(b, 2) - (4 * a * c);
    float t_max = 0.0;
    if (determinant > 0) {
      float[] roots = {( -b + sqrt(determinant) ) / (2 * a), ( -b - sqrt(determinant) ) / (2 * a)};
      
      if (roots[0] > 0 && roots[1] > 0) {
        t_max = roots[0] < roots[1] ? roots[0] : roots[1];
      } else if (roots[0] < 0 && roots[1] > 0) {
         t_max = roots[1];
      } else if (roots[0] > 0 && roots[1] < 0) {
        t_max = roots[0];
      } else { 
        t_max = 0.0;
      }

      if (t_max > 0.0 && t_max < minRoot) {
        minRoot = t_max;
        
        PVector intersect = PVector.add(ray.origin, PVector.mult(ray.direction, minRoot));
          
        PVector hitNormal = PVector.sub(intersect, ellipsoid.pos);
        hitNormal = new PVector(hitNormal.x / pow(ellipsoid.radi.x, 2), 
          hitNormal.y / pow(ellipsoid.radi.y, 2),
          hitNormal.z / pow(ellipsoid.radi.z, 2)); 
        hitNormal = hitNormal.normalize();
        hit = new Hit(ellipsoid, minRoot, intersect, hitNormal, ray);
      }
    }
  }
  return hit;
}

// lights
class Light { 
  PVector pos, colors;

  Light(float x, float y, float z, float r, float g, float b) {  
    pos = new PVector(x, y, z);
    colors = new PVector(r, g, b);
  }
}

// surfaces
class Surface {  
  PVector diffuse;
  PVector ambient;
  PVector specular;
  float phongExp, Krefl;

  Surface(float Cdr, float Cdg, float Cdb, float Car, float Cag, float Cab, float Csr, float Csg, float Csb, float P, float K) {
    diffuse = new PVector(Cdr, Cdg, Cdb);
    ambient = new PVector(Car, Cag, Cab);
    specular = new PVector(Csr, Csg, Csb);

    phongExp = P;
    Krefl = K;
  }
}

// shapes
abstract class Shape {
  PVector pos, diffuse, ambient, specular;
  float phongExp, Krefl;
}

class Disk extends Shape {
  PVector normal;
  float radius;

  Disk(float x, float y, float z, float nx, float ny, float nz, float r, Surface surface) {  
    pos = new PVector(x, y, z);
    diffuse = surface.diffuse;
    ambient = surface.ambient;
    specular = surface.specular;
    phongExp = surface.phongExp;
    Krefl = surface.Krefl;

    normal = new PVector(nx, ny, nz);
    radius = r;
  }
}

class Ellipsoid extends Shape {
  PVector radi;

  float x, y, z, rx, ry, rz;
  Ellipsoid (float x, float y, float z, float rx, float ry, float rz, Surface surface) {  
    pos = new PVector(x, y, z);
    diffuse = surface.diffuse;
    ambient = surface.ambient;
    specular = surface.specular;
    phongExp = surface.phongExp;
    Krefl = surface.Krefl;

    radi = new PVector(rx, ry, rz);
  }
}
