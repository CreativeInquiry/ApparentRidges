// Processing Example for ApparentRidges, the 3D edge detection toolkit
// INSTRUCTIONS

// A pre-built Processing library can be found at processing/ApparentRidges
// just place it in your Processing libraries folder.
// for example,
// cp -r processing/ApparentRidges ~/Documents/Processing/libraries
// done!

// Alternatively, you can also download Eclipse, set up from Processing
// library template, copy paste the java source code from java/ApparentRidges,
// and ... compile. But you'll basically get the same .jar file, so why bother!

// ===========================================================================

// import the library
import apparentridges.*;

Render render;
PImage pg;

void setup(){
  size(800,800);
  Mesh mesh = OBJParser.fromFile(sketchPath("../../testdata/dragon_100K.obj"));
  render = new Render(mesh,width,height);
  render.focal = 500;
  render.transform(Util.matTrsl(0,-0.1,0));
  render.transform(Util.matScal(100,100,100));
  render.transform(Util.matTrsl(0,0,20));
  
  render.mesh.precompute(true,true);
  
  render.didPrecompute=true; // tell renderer to not precompute the mesh again
                             // when trying to calculate ridges, since we just
                             // called precompute explicitly.
  
  // retrieve a curvature map too (not required, just for showing off)
  double[] dat = PixelMap.curvature(render);
  pg = createImage(width,height,RGB);
  pg.loadPixels();
  for (int i = 0; i < dat.length; i+=2){
    pg.pixels[i/2] = color((int)(dat[i]*255+128),(int)(dat[i+1]*255+128),128);
  }
  pg.updatePixels();
}
void draw(){
  background(255);
  if (keyPressed){
    image(pg,0,0);
  }
  noFill();
  stroke(0);
 
  render.clear();
  render.apparentRidges(1-(float)mouseY/(float)height,2);
  for (int i = 0; i < render.lines.length; i++){
    Line l = render.lines.__get(i);
    line((float)l.x1,(float)l.y1,(float)l.x2,(float)l.y2);
  }
  
  fill(0);
  noStroke();
  text("use mouseY to control threshold, hold key to show curvature map, click to save svg",5,12);
}


void mousePressed(){
  render.apparentRidges((float)mouseY/(float)height,2);
  render.buildPolylines(null);
  String out = SVGWriter.polylines(render,null);
  Util.writeFile(sketchPath("output.svg"),out);
  println("saved!");
  
}
