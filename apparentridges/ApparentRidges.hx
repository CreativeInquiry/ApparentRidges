/***************************************************
 * Apparent Ridges 
 *
 * -- Render line drawings of 3D meshes.
 *
 * Main Algorithm based on RTSC (C++): https://rtsc.cs.princeton.edu
 * Apparent Ridges Paper: http://people.csail.mit.edu/tjudd/apparentLines.pdf
 *
 * This haxe port by:
 * Lingdong Huang 2020, License: GPL v2
 * Made for STUDIO for Creative Inquiry at CMU
 *
 * haxe --version
 * 4.1.3
 * See Makefile for build instructions
 ***************************************************/

package apparentridges;

import haxe.ds.Vector;

@:expose
class Util {
  public static inline function rndInt(x : Float) : Int{
    return Std.int(Math.round(x));
  }
  public static inline function min(a:Int, b:Int) : Int{
    return (a>b)?b:a;
  }

  public static inline function sq(x : Float) : Float{
    return x * x;
  }
  // i + 1 and i - 1 modulo 3
  // This way of computing it tends to be faster than using %
  public static inline function nextMod3(i : Int) : Int{
    return ((i) < 2 ? (i) + 1 : (i) - 2);
  }
  public static inline function prevMod3(i : Int) : Int{
    return ((i) > 0 ? (i) - 1 : (i) + 2);
  }
  // Area-weighted triangle face normal
  public static inline function trinorm(v0 : Vec3, v1 : Vec3, v2 : Vec3) : Vec3{
    return Vec3.cross((v1 - v0), (v2 - v0)) * 0.5;
  }

  // LDL^T decomposition of a symmetric positive definite matrix (and some
  // other symmetric matrices, but fragile since we don't do pivoting).
  // Like Cholesky, but no square roots, which is important for small N.
  // Reads diagonal and upper triangle of matrix A.
  // On output, lower triangle of A holds LD, while rdiag holds D^-1.
  // Algorithm from Golub and van Loan.
  public static function ldltdc(A : Array<Array<Float>>, rdiag : Array<Float>) : Bool{
    // Special case for small N
    var N : Int = rdiag.length;
    if (N < 1) {
      return false;
    } else if (N <= 3) {
      var d0 : Float = A[0][0];
      rdiag[0] = 1 / d0;
      if (N == 1)
        return (d0 != 0);
      A[1][0] = A[0][1];
      var l10 : Float = rdiag[0] * A[1][0];
      var d1 : Float = A[1][1] - l10 * A[1][0];
      rdiag[1] = 1 / d1;
      if (N == 2)
        return (d0 != 0 && d1 != 0);
      var d2 : Float = A[2][2] - rdiag[0] * sq(A[2][0]) - rdiag[1] * sq(A[2][1]);
      rdiag[2] = 1 / d2;
      A[2][0] = A[0][2];
      A[2][1] = A[1][2] - l10 * A[2][0];
      return (d0 != 0 && d1 != 0 && d2 != 0);
    }
    var v : Array<Float> = []; //[N-1]
    for (i in 0...N) {
      for (k in 0...i)
        v[k] = A[i][k] * rdiag[k];
      for (j in i...N) {
        var sum : Float = A[i][j];
        for (k in 0...i)
          sum -= v[k] * A[j][k];
        if (i == j) {
          if (sum == 0)
            return false;
          rdiag[i] = 1 / sum;
        } else {
          A[j][i] = sum;
        }
      }
    }

    return true;
  }
  // Solve Ax=b after ldltdc.  x is allowed to be the same as b.
  public static function ldltsl(A : Array<Array<Float>>, 
                                rdiag : Array<Float>,
                                b : Array<Float>,
                                x : Array<Float>){
    var N : Int = rdiag.length;
    for (i in 0...N) {
      var sum : Float = b[i];
      for (k in 0...i)
        sum -= A[i][k] * x[k];
      x[i] = sum * rdiag[i];
    }
    for (_i in 0...N) {
      var i : Int = N-_i-1;
      var sum : Float = 0;
      for (k in (i+1)...N)
        sum += A[k][i] * x[k];
      x[i] -= sum * rdiag[i];
    }
  }

  // Compute largest eigenvalue and associated eigenvector of a
  // symmetric 2x2 matrix.  Solves characteristic equation.
  // Inputs: three elements of matrix (upper-left, diag, lower-right)
  // Outputs: largest (in magnitude) eigenvector/value
  public static function largestEig2x2(m1:Float,m12:Float,m2:Float) : Vec3{
    var l1 : Float = 0.5 * (m1+m2);
    
    // The result of the below sqrt is positive, so to get the largest
    // eigenvalue we add it if we were positive already, else subtract
    if (l1 > 0.0)
      l1 += Math.sqrt(Util.sq(m12) + 0.25 * Util.sq(m2-m1));
    else
      l1 -= Math.sqrt(Util.sq(m12) + 0.25 * Util.sq(m2-m1));

    // Find corresponding eigenvector
    var e1 = new Vec3(m2 - l1, -m12, 0);
    e1.normalize();
    e1.z = l1;
    return e1;
  }

  public static inline function matIden() : Array<Float> {return [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];}
  public static inline function matRotx(a:Float): Array<Float> {return [1,0,0,0, 0,Math.cos(a),-Math.sin(a),0, 0,Math.sin(a),Math.cos(a),0, 0,0,0,1];}
  public static inline function matRoty(a:Float): Array<Float> {return [Math.cos(a),0,Math.sin(a),0, 0,1,0,0, -Math.sin(a),0,Math.cos(a),0, 0,0,0,1];}
  public static inline function matRotz(a:Float): Array<Float> {return [Math.cos(a),-Math.sin(a),0,0, Math.sin(a),Math.cos(a),0,0, 0,0,1,0, 0,0,0,1];}
  public static inline function matTrsl(x:Float,y:Float,z:Float) : Array<Float> {return [1,0,0,x, 0,1,0,y, 0,0,1,z, 0,0,0,1];}
  public static inline function matScal(x:Float,y:Float,z:Float) : Array<Float> {return [x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1];}
  public static inline function matMult(A:Array<Float>,B:Array<Float>) : Array<Float> {
    return [A[0]*B[0]+A[1]*B[4]+A[2]*B[8]+A[3]*B[12],A[0]*B[1]+A[1]*B[5]+A[2]*B[9]+A[3]*B[13],A[0]*B[2]+A[1]*B[6]+A[2]*B[10]+A[3]*B[14],
      A[0]*B[3]+A[1]*B[7]+A[2]*B[11]+A[3]*B[15],A[4]*B[0]+A[5]*B[4]+A[6]*B[8]+A[7]*B[12],A[4]*B[1]+A[5]*B[5]+A[6]*B[9]+A[7]*B[13],
      A[4]*B[2]+A[5]*B[6]+A[6]*B[10]+A[7]*B[14],A[4]*B[3]+A[5]*B[7]+A[6]*B[11]+A[7]*B[15],A[8]*B[0]+A[9]*B[4]+A[10]*B[8]+A[11]*B[12],
      A[8]*B[1]+A[9]*B[5]+A[10]*B[9]+A[11]*B[13],A[8]*B[2]+A[9]*B[6]+A[10]*B[10]+A[11]*B[14],A[8]*B[3]+A[9]*B[7]+A[10]*B[11]+A[11]*B[15],
      A[12]*B[0]+A[13]*B[4]+A[14]*B[8]+A[15]*B[12],A[12]*B[1]+A[13]*B[5]+A[14]*B[9]+A[15]*B[13],A[12]*B[2]+A[13]*B[6]+A[14]*B[10]+A[15]*B[14],
      A[12]*B[3]+A[13]*B[7]+A[14]*B[11]+A[15]*B[15]];
  }
  public static inline function matTrfm(A:Array<Float>,v:Vec3) : Vec3 {
    return new Vec3((A[0]*v[0]+A[1]*v[1]+A[2]*v[2]+A[3])/(A[12]*v[0]+A[13]*v[1]+A[14]*v[2]+A[15]),
      (A[4]*v[0]+A[5]*v[1]+A[6]*v[2]+A[7])/(A[12]*v[0]+A[13]*v[1]+A[14]*v[2]+A[15]),(A[8]*v[0]+A[9]*v[1]+A[10]*v[2]+A[11])/(A[12]*v[0]+A[13]*v[1]+A[14]*v[2]+A[15])
    );}
  public static inline function matProj(f:Float,v:Vec3) : Vec3 {return new Vec3(f*v[0]/v[2],f*v[1]/v[2],0);}
  
  public static inline function uniformHemisphereSampler () : Vec3{
    var Xi1 : Float = Math.random();
    var Xi2 : Float = Math.random();
    var theta : Float = Math.acos(Xi1);
    var phi : Float = 2*Math.PI*Xi2;
    var xs : Float = Math.sin(theta)*Math.cos(phi);
    var ys : Float = Math.sin(theta)*Math.sin(phi);
    var zs : Float = Math.cos(theta);
    var v : Vec3 = new Vec3(xs,ys,zs);
    return v;
  }

#if sys
  public static function writeFile(filename : String, content : String){
    sys.io.File.saveContent(filename,content);
  }
#end
}

@:expose
abstract Vec3 (haxe.ds.Vector<Float>) {
  public var x(get, set):Float;
  public var y(get, set):Float;
  public var z(get, set):Float;

  public inline function new(_x:Float,_y:Float,_z:Float){
    this = new haxe.ds.Vector<Float>(3);
    this[0] = _x;
    this[1] = _y;
    this[2] = _z;
  }
  private inline function get_x():Float {  return this[0]; }
  private inline function get_y():Float {  return this[1]; }
  private inline function get_z():Float {  return this[2]; }
  private inline function set_x(v:Float):Float {  return this[0]=v; }
  private inline function set_y(v:Float):Float {  return this[1]=v; }
  private inline function set_z(v:Float):Float {  return this[2]=v; }

  public inline function copy() : Vec3{
    return new Vec3(x,y,z);
  }
  public inline function assign(v : Vec3){
    x = v.x; y = v.y; z = v.z;
  }

  public static inline function cross(v1:Vec3, v2:Vec3) : Vec3 {
    return new Vec3(
      v1.y*v2.z-v1.z*v2.y,
      v1.z*v2.x-v1.x*v2.z,
      v1.x*v2.y-v1.y*v2.x
    );
  }
  public static inline function dot(v1:Vec3, v2:Vec3) : Float {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  }
  public static inline function dist2(v1:Vec3, v2:Vec3) : Float{
    return Util.sq(v1.x-v2.x)+Util.sq(v1.y-v2.y)+Util.sq(v1.z-v2.z);
  }
  public static inline function dist(v1:Vec3,v2:Vec3):Float{
    return Math.sqrt(Vec3.dist2(v1,v2));
  }

  public inline function len():Float {
    return Math.sqrt(x*x+y*y+z*z);
  }
  public inline function len2():Float{
    return x*x+y*y+z*z;
  }

  public inline function normalize(){
    var l : Float = len();
    if (l > 0){
      l = 1/l;
      x *= l;
      y *= l;
      z *= l;
    }else{
      x = 0;
      y = 0;
      z = 1;
    }
  }
  public inline function scale(s:Float){
    x *= s;
    y *= s;
    z *= s;
  }

  @:arrayAccess
  public inline function get(i:Int) {
    return this[i];
  }
  @:arrayAccess
  public inline function set(i:Int,v:Float):Float {
    this[i] = v;
    return v;
  }

  @:op(A+B)
  public inline function add(rhs:Vec3) : Vec3 {
    return new Vec3(x+rhs.x,y+rhs.y,z+rhs.z);
  }
  @:op(A-B)
  public inline function sub(rhs:Vec3) : Vec3 {
    return new Vec3(x-rhs.x,y-rhs.y,z-rhs.z);
  }
  @:op(A*B)
  public inline function mul(rhs:Vec3) : Vec3 {
    return new Vec3(x*rhs.x,y*rhs.y,z*rhs.z);
  }
  @:op(A*B)
  public inline function mulf(rhs:Float) : Vec3 {
    return new Vec3(x*rhs,y*rhs,z*rhs);
  }

}

@:expose
abstract Face (haxe.ds.Vector<Int>) {
  public inline function new(a:Int, b:Int, c:Int){
    this = new haxe.ds.Vector<Int>(3);
    this[0]=a;
    this[1]=b;
    this[2]=c;
  }
  @:arrayAccess
  public inline function get(i:Int) {
    return this[i];
  }
  @:arrayAccess
  public inline function set(i:Int,v:Int):Int {
    this[i] = v;
    return v;
  }
  public inline function indexOf(v:Int):Int{
    if (this[0]==v){return 0;}
    if (this[1]==v){return 1;}
    if (this[2]==v){return 2;}
    return -1;
  }
}

@:expose
class Ridge {
  public var A : Vec3;
  public var B : Vec3;
  public var strengthA : Float;
  public var strengthB : Float;
  public function new(a:Vec3,sa:Float,b:Vec3,sb:Float){
    A = a; strengthA = sa;
    B = b; strengthB = sb;
  }
}

@:expose
class BSphere{
  public var r : Float;
  public var o : Vec3;
  public function new(){}
}

@:expose
class Mesh {
  public var vertices : Array<Vec3> = [];
  public var faces : Array<Face> = [];
  public var normals : Array<Vec3> = [];
  public var curv1 : Array<Float> = [];
  public var curv2 : Array<Float> = [];
  public var pdir1 : Array<Vec3> = [];
  public var pdir2 : Array<Vec3> = [];
  public var pointAreas : Array<Float> = [];
  public var cornerAreas : Array<Vec3> = [];
  public var adjacentFaces : Array<Array<Int>> = [];

  public var ndotv : Array<Float> = [];
  public var t1q1 : Array<Vec3> = [];
  public var Dt1q1 : Array<Float> = [];

  public var bsphere : BSphere;
  public var featureSize : Float;

  public var bvh : BVHTree;

  public function new(){

  }

  public function precompute(doBVH:Bool=true,verb:Bool=false){
    if(verb)trace("computing normals...");
    computeNormals();
    if(verb)trace("computing point areas...");
    computePointAreas();
    if(verb)trace("computing adjacent faces...");
    computeAdjacentFaces();
    if(verb)trace("computing curvatures...");
    computeCurvatures();
    if(verb)trace("computing bounding sphere...");
    computeBSphere();
    if(verb)trace("computing feature size...");
    computeFeatureSize();
    if(verb)trace("computing bounding volume hierarchy...");
    if (doBVH){computeBVH();}else{computeBVHTrivial();}
    if(verb)trace("pre-computation finished.");
  }

  public function computeNormals(){
    normals = [for (i in 0...vertices.length) new Vec3(0,0,0)];
    for (f in faces){
      var p0 : Vec3 = vertices[f[0]];
      var p1 : Vec3 = vertices[f[1]];
      var p2 : Vec3 = vertices[f[2]];
      var a  : Vec3 = p0-p1;
      var b  : Vec3 = p1-p2;
      var fn : Vec3 = Vec3.cross(a,b);
      normals[f[0]] += fn;
      normals[f[1]] += fn;
      normals[f[2]] += fn;
    }
    for (n in normals){
      n.normalize();
    }
  }

  public function computeBSphere(){
    bsphere = new BSphere();

    function farthestVertexAlong(dir : Vec3){
      var nv : Int = vertices.length;
      var far : Int = 0;
      var far_dot = Vec3.dot(vertices[0],dir);
      for (i in 1...nv){
        var my_dot : Float = Vec3.dot(vertices[i],dir);
        if (my_dot > far_dot){
          far = i;
          far_dot = my_dot;
        }
      }
      return far;
    }
    var best_min : Vec3 = new Vec3(0,0,0);
    var best_max : Vec3 = new Vec3(0,0,0);
    var dirs : Array<Vec3> = [
      new Vec3(1,0,0),
      new Vec3(0,1,0),
      new Vec3(0,0,1),
      new Vec3(1,1,1),
      new Vec3(1,-1,1),
      new Vec3(1,-1,-1),
      new Vec3(1,1,1),
    ];
    for (d in dirs){
      var p1 : Vec3 = vertices[farthestVertexAlong(d*(-1.0))];
      var p2 : Vec3 = vertices[farthestVertexAlong(d)];
      if (Vec3.dist2(p1,p2)>Vec3.dist2(best_min,best_max)){
        best_min = p1;
        best_max = p2;
      }
    }
    bsphere.o = (best_min + best_max) * 0.5;
    bsphere.r = Vec3.dist(best_min, best_max) * 0.5;
    var r2 = Util.sq(bsphere.r);
    // Expand bsphere to contain all points
    for (i in 0...vertices.length){
      var d2 : Float = Vec3.dist2(vertices[i],bsphere.o);
      if (d2 <= r2){
        continue;
      }
      var d : Float = Math.sqrt(d2);
      bsphere.r = 0.5 * (bsphere.r + d);
      r2 = Util.sq(bsphere.r);
      bsphere.o-=vertices[i];
      bsphere.o*=bsphere.r*(1.0/d);
      bsphere.o+=vertices[i];
    }

  }

  // Compute a "feature size" for the mesh: computed as 1% of
  // the reciprocal of the 10-th percentile curvature
  public function computeFeatureSize(){
    var nv : Int = curv1.length;
    var nsamp = Util.min(nv, 500);
    var samples : Array<Float> = [];
    var s : Int = 79;//bbs
    var p : Int = 103;
    var q : Int = 211;
    var m : Int = p*q;
    for (i in 0...nsamp){
      var ind : Int = Std.int(nv*(s/m));
      s = (s*s) % m;
      samples.push(Math.abs(curv1[ind]));
      samples.push(Math.abs(curv2[ind]));
    }
    var frac : Float = 0.1;
    var mult : Float = 0.01;
    var max_feat_size = 0.05 * bsphere.r;
    var which : Int = Std.int(frac * samples.length);
    samples.sort(function(a : Float, b : Float) : Int {
      if(a < b) return -1;
      else if(a > b) return 1;
      else return 0;
    }); // todo: std::nth_element
    featureSize = Math.min(mult/samples[which],max_feat_size);
  }

  public function computeAdjacentFaces(){
    adjacentFaces = [for (i in 0...vertices.length) new Array<Int>()];
    for (i in 0...faces.length){
      for (j in 0...3){
        adjacentFaces[faces[i][j]].push(i);
      }
    }
  }

  public function getFaceEdges(f:Face) : haxe.ds.Vector<Vec3>{
      var e : haxe.ds.Vector<Vec3> = new haxe.ds.Vector<Vec3>(3);
      e[0] = vertices[f[2]]-vertices[f[1]];
      e[1] = vertices[f[0]]-vertices[f[2]];
      e[2] = vertices[f[1]]-vertices[f[0]];
      return e;
  }

  public function computePointAreas(){
    var nf : Int = faces.length;
    var nv : Int = vertices.length;

    pointAreas = [for (i in 0...nv) 0];
    cornerAreas = [for (i in 0...nf) new Vec3(0,0,0)];

    for (i in 0...nf){
      var e : haxe.ds.Vector<Vec3> = getFaceEdges(faces[i]);
      // Compute corner weights
      var area : Float = 0.5 * Vec3.cross(e[0],e[1]).len();
      var l2 : Array<Float> = [e[0].len2(), e[1].len2(), e[2].len2()];
      // Barycentric weights of circumcenter
      var bcw : Array<Float> = [
        l2[0] * (l2[1] + l2[2] - l2[0]),
        l2[1] * (l2[2] + l2[0] - l2[1]),
        l2[2] * (l2[0] + l2[1] - l2[2]) 
      ];
      if (bcw[0] <= 0) {
        cornerAreas[i][1] = -0.25 * l2[2] * area / Vec3.dot(e[0],e[2]);
        cornerAreas[i][2] = -0.25 * l2[1] * area / Vec3.dot(e[0],e[1]);
        cornerAreas[i][0] = area - cornerAreas[i][1] - cornerAreas[i][2];
      } else if (bcw[1] <= 0.0) {
        cornerAreas[i][2] = -0.25 * l2[0] * area / Vec3.dot(e[1],e[0]);
        cornerAreas[i][0] = -0.25 * l2[2] * area / Vec3.dot(e[1],e[2]);
        cornerAreas[i][1] = area - cornerAreas[i][2] - cornerAreas[i][0];
      } else if (bcw[2] <= 0.0) {
        cornerAreas[i][0] = -0.25 * l2[1] * area / Vec3.dot(e[2],e[1]);
        cornerAreas[i][1] = -0.25 * l2[0] * area / Vec3.dot(e[2],e[0]);
        cornerAreas[i][2] = area - cornerAreas[i][0] - cornerAreas[i][1];
      } else {
        var scale : Float = 0.5 * area / (bcw[0] + bcw[1] + bcw[2]);
        for (j in 0...3)
          cornerAreas[i][j] = scale * (bcw[Util.nextMod3(j)] +
                                       bcw[Util.prevMod3(j)]);
      }
      pointAreas[faces[i][0]] += cornerAreas[i][0];
      pointAreas[faces[i][1]] += cornerAreas[i][1];
      pointAreas[faces[i][2]] += cornerAreas[i][2];
    }
  }

  // Rotate a coordinate system to be perpendicular to the given normal
  private static function rotCoordSys(old_u : Vec3, old_v : Vec3, new_norm : Vec3,
                                      new_u : Vec3, new_v : Vec3){
    new_u.assign(old_u);
    new_v.assign(old_v);
    var old_norm : Vec3 = Vec3.cross(old_u,old_v);
    var ndot : Float = Vec3.dot(old_norm,new_norm);
    if ((ndot <= -1)) {
      new_u.scale(-1);
      new_v.scale(-1);
      return;
    }
    // Perpendicular to old_norm and in the plane of old_norm and new_norm
    var perp_old : Vec3 = new_norm - old_norm*ndot;

    // Perpendicular to new_norm and in the plane of old_norm and new_norm
    // vec perp_new = ndot * new_norm - old_norm;

    // perp_old - perp_new, with normalization constants folded in
    var dperp : Vec3 = (old_norm + new_norm) * (1 / (1 + ndot));

    // Subtracts component along perp_old, and adds the same amount along
    // perp_new.  Leaves unchanged the component perpendicular to the
    // plane containing old_norm and new_norm.
    new_u -= dperp * Vec3.dot(new_u,perp_old);
    new_v -= dperp * Vec3.dot(new_v,perp_old);
  }

  // Reproject a curvature tensor from the basis spanned by old_u and old_v
  // (which are assumed to be unit-length and perpendicular) to the
  // new_u, new_v basis.
  private static function projCurv(old_u : Vec3, old_v : Vec3,
                                   old_ku : Float, old_kuv : Float, old_kv : Float,
                                   new_u : Vec3, new_v : Vec3) : Vec3{
    var r_new_u : Vec3 = new Vec3(0,0,0);
    var r_new_v : Vec3 = new Vec3(0,0,0);
    rotCoordSys(new_u,new_v,Vec3.cross(old_u,old_v),r_new_u,r_new_v);
    var u1 : Float = Vec3.dot(r_new_u,old_u);
    var v1 : Float = Vec3.dot(r_new_u,old_v);
    var u2 : Float = Vec3.dot(r_new_v,old_u);
    var v2 : Float = Vec3.dot(r_new_v,old_v);
    var new_ku  : Float = old_ku * u1*u1 + old_kuv * (2     * u1*v1) + old_kv * v1*v1;
    var new_kuv : Float = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
    var new_kv  : Float = old_ku * u2*u2 + old_kuv * (2     * u2*v2) + old_kv * v2*v2;
    return new Vec3(new_ku,new_kuv,new_kv);
  }

  // Given a curvature tensor, find principal directions and curvatures
  // Makes sure that pdir1 and pdir2 are perpendicular to normal
  private static function diagonalizeCurv(old_u : Vec3, old_v : Vec3,
                                          ku : Float, kuv : Float, kv : Float,
                                          new_norm : Vec3, 
                                          pd1 : Vec3, pd2 : Vec3, k1k2 : Vec3){
    var r_old_u : Vec3 = new Vec3(0,0,0); 
    var r_old_v : Vec3 = new Vec3(0,0,0); 
    rotCoordSys(old_u, old_v, new_norm, r_old_u, r_old_v);

    var c : Float = 1; var s : Float = 0; var tt : Float = 0;
    if (kuv != 0.0) {
      // Jacobi rotation to diagonalize
      var h : Float = 0.5 * (kv - ku) / kuv;
      tt = (h < 0.0) ?
        1.0 / (h - Math.sqrt(1.0 + h*h)) :
        1.0 / (h + Math.sqrt(1.0 + h*h));
      c = 1.0 / Math.sqrt(1.0 + tt*tt);
      s = tt * c;
    }

    var k1 : Float = ku - tt * kuv;
    var k2 : Float = kv + tt * kuv;

    if (Math.abs(k1) >= Math.abs(k2)) {
      k1k2.x = k1;
      k1k2.y = k2;
      pd1.assign(r_old_u*c - r_old_v*s);
    } else {
      k1k2.x = k2;
      k1k2.y = k1;
      pd1.assign(r_old_u*s + r_old_v*c);
    }
    pd2.assign(Vec3.cross(new_norm,pd1));
  }


  // compute principal curvatures and directions
  public function computeCurvatures(){

    var nv : Int = vertices.length;
    var nf : Int = faces.length;

    curv1 = [for (i in 0...nv) 0];
    curv2 = [for (i in 0...nv) 0];
    pdir1 = [for (i in 0...nv) new Vec3(0,0,0)];
    pdir2 = [for (i in 0...nv) new Vec3(0,0,0)];
    var curv12 : Array<Float> = [for (i in 0...nv) 0];

    // Set up an initial coordinate system per vertex
    for (i in 0...nf){
      pdir1[faces[i][0]] = vertices[faces[i][1]] -
                           vertices[faces[i][0]];
      pdir1[faces[i][1]] = vertices[faces[i][2]] -
                           vertices[faces[i][1]];
      pdir1[faces[i][2]] = vertices[faces[i][0]] -
                           vertices[faces[i][2]];
    }
    for (i in 0...nv){
      pdir1[i] = Vec3.cross(pdir1[i],normals[i]);
      pdir1[i].normalize();
      pdir2[i] = Vec3.cross(normals[i],pdir1[i]);
    }

    // Compute curvature per-face
    for (i in 0...nf){
      var f : Face = faces[i];
      var e : haxe.ds.Vector<Vec3> = getFaceEdges(f);
      
      // N-T-B coordinate system per face
      var t : Vec3 = e[0].copy();
      t.normalize();
      var n : Vec3 = Vec3.cross(e[0],e[1]);
      var b : Vec3 = Vec3.cross(n,t);
      b.normalize();
      // Estimate curvature based on variation of normals
      // along edges
      var m : Array<Float> = [0,0,0];
      var w : Array<Array<Float>> = [[0,0,0],[0,0,0],[0,0,0]];
      for (j in 0...3){
        var u : Float = Vec3.dot(e[j],t);
        var v : Float = Vec3.dot(e[j],b);
        w[0][0] += u*u;
        w[0][1] += u*v;
        w[2][2] += v*v;
        // The below are computed once at the end of the loop
        // w[1][1] += v*v + u*u;
        // w[1][2] += u*v;
        var dn : Vec3 = normals[f[Util.prevMod3(j)]] -
                        normals[f[Util.nextMod3(j)]];
        var dnu : Float = Vec3.dot(dn,t);
        var dnv : Float = Vec3.dot(dn,b);
        m[0] += dnu*u;
        m[1] += dnu*v + dnv*u;
        m[2] += dnv*v;
      }
      w[1][1] = w[0][0] + w[2][2];
      w[1][2] = w[0][1];
      // Least squares solution
      var diag : Array<Float> = [0,0,0];
      if (!Util.ldltdc(w, diag)) {
        continue;
      }
      Util.ldltsl(w,diag,m,m);
      // Push it back out to the vertices
      for (j in 0...3){
        var vj : Int = f[j];
        var ccc : Vec3 = projCurv(t,b,m[0],m[1],m[2],pdir1[vj],pdir2[vj]);
        var c1 : Float= ccc.x;
        var c12: Float= ccc.y;
        var c2 : Float= ccc.z;
        var wt : Float = cornerAreas[i][j]/pointAreas[vj];
        curv1[vj] += wt * c1;
        curv12[vj] += wt * c12;
        curv2[vj] += wt * c2;

      }
    }
    for (i in 0...nv){
      var c1c2 : Vec3 = new Vec3(0,0,0);
      diagonalizeCurv(
        pdir1[i],pdir2[i],
        curv1[i],curv12[i],curv2[i],
        normals[i],pdir1[i],pdir2[i],
        c1c2
      );
      curv1[i]=c1c2.x;
      curv2[i]=c1c2.y;
    }

  }


  private function computeVertViewDepCurv(i:Int,ndotv:Float,u2:Float,
                                          uv:Float,v2:Float) : Vec3{
    // Find the entries in Q = S * P^-1
    //                       = S + (sec theta - 1) * S * w * w^T
    var sectheta_minus1 : Float = 1.0 / Math.abs(ndotv) - 1.0;
    var Q11: Float = curv1[i] * (1.0 + sectheta_minus1 * u2);
    var Q12: Float = curv1[i] * (      sectheta_minus1 * uv);
    var Q21: Float = curv2[i] * (      sectheta_minus1 * uv);
    var Q22: Float = curv2[i] * (1.0 + sectheta_minus1 * v2);

    // Find the three entries in the (symmetric) matrix Q^T Q
    var QTQ1 : Float = Q11 * Q11 + Q21 * Q21;
    var QTQ12: Float = Q11 * Q12 + Q21 * Q22;
    var QTQ2 : Float = Q12 * Q12 + Q22 * Q22;

    // Compute eigenstuff
    return Util.largestEig2x2(QTQ1, QTQ12, QTQ2);
  }

  // Compute D_{t_1} q_1 - the derivative of max view-dependent curvature
  // in the principal max view-dependent curvature direction.
  private function computeVertDt1q1(i : Int, ndotv : Float,
                                    t1q1 : Array<Vec3>) : Float {
    var v0 : Vec3 = vertices[i];
    var this_viewdep_curv : Float = t1q1[i].z;
    var world_t1 : Vec3 =  pdir1[i] * t1q1[i][0] + pdir2[i] * t1q1[i][1];
    var world_t2 : Vec3 = Vec3.cross(normals[i], world_t1);
    var v0_dot_t2 : Float = Vec3.dot(v0,world_t2);

    var Dt1q1 : Float = 0.0;
    var n : Int= 0;

    var naf : Int = adjacentFaces[i].length;

    for (j in 0...naf) {
      // We're in a triangle adjacent to the vertex of interest.
      // The current vertex is v0 - let v1 and v2 be the other two
      var f  : Int=  adjacentFaces[i][j];
      var ind: Int = faces[f].indexOf(i);
      var i1 : Int = faces[f][Util.nextMod3(ind)];
      var i2 : Int = faces[f][Util.prevMod3(ind)];
      var v1 : Vec3 = vertices[i1];
      var v2 : Vec3 = vertices[i2];

      // Find the point p on the segment between v1 and v2 such that
      // its vector from v0 is along t1, i.e. perpendicular to t2.
      // Linear combination: p = w1*v1 + w2*v2, where w2 = 1-w1
      var v1_dot_t2 : Float = Vec3.dot(v1 , world_t2);
      var v2_dot_t2 : Float = Vec3.dot(v2 , world_t2);
      var w1 : Float = (v2_dot_t2 - v0_dot_t2) / (v2_dot_t2 - v1_dot_t2);

      // If w1 is not in [0..1) then we're not interested.  
      // Incidentally, the computation of w1 can result in infinity,
      // but the comparison should do the right thing...
      if (w1 < 0.0 || w1 >= 1.0)
        continue;
      
      // Construct the opposite point
      var w2 : Float = 1.0 - w1;
      var p : Vec3 = v1 * w1 + v2 * w2;

      // And interpolate to find the view-dependent curvature at that point
      var interp_viewdep_curv : Float = w1 * t1q1[i1].z + w2 * t1q1[i2].z;

      // Finally, take the *projected* view-dependent curvature derivative
      var proj_dist : Float = Vec3.dot((p - v0),world_t1);
      proj_dist *= Math.abs(ndotv);
      Dt1q1 += (interp_viewdep_curv - this_viewdep_curv) / proj_dist;
      n++;

      // To save time, quit as soon as we have two estimates
      // (that's all we're going to get, anyway)
      if (n == 2) {
        Dt1q1 *= 0.5;
        return Dt1q1;
      }
    }
    return Dt1q1;
  }

  // Draw part of a ridge/valley curve on one triangle face.  v0,v1,v2
  // are the indices of the 3 vertices; this function assumes that the
  // curve connects points on the edges v0-v1 and v1-v2
  // (or connects point on v0-v1 to center if to_center is true)
  public function segmentApparentRidge(
      v0 : Int, v1 : Int, v2 : Int,
      emax0 : Float, emax1 : Float, emax2 : Float,
      kmax0 : Float, kmax1 : Float, kmax2 : Float,
      tmax0 : Vec3, tmax1 : Vec3, tmax2 : Vec3,
      thresh : Float, to_center : Bool, do_test : Bool
    ) : Ridge {

    // Interpolate to find ridge/valley line segment endpoints
    // in this triangle and the curvatures there
    var w10 : Float = Math.abs(emax0) / (Math.abs(emax0) + Math.abs(emax1));
    var w01 : Float = 1.0 - w10;
    var p01 : Vec3 = vertices[v0] * w01 + vertices[v1] * w10;
    var k01 : Float = Math.abs(w01 * kmax0 + w10 * kmax1);

    var p12 : Vec3;
    var k12 : Float;
    if (to_center){
      // Connect first point to center of triangle
      p12 =(vertices[v0] +
            vertices[v1] +
            vertices[v2]) * (1.0/3.0);
      k12 = Math.abs(kmax0 + kmax1 + kmax2) / 3.0;
    }else{
      // Connect first point to second one (on next edge)
      var w21 : Float = Math.abs(emax1) / (Math.abs(emax1) + Math.abs(emax2));
      var w12 : Float = 1.0 - w21;
      p12 = vertices[v1] * w12 + vertices[v2] * w21;
      k12 = Math.abs(w12 * kmax1 + w21 * kmax2);
    }
    // Don't draw below threshold
    k01 -= thresh;
    if (k01 < 0.0)
      k01 = 0.0;
    k12 -= thresh;
    if (k12 < 0.0)
      k12 = 0.0;

    // Skip lines that you can't see...
    if (k01 == 0.0 && k12 == 0.0)
      return null;

    // Perform test: do the tmax-es point *towards* the segment? (Fig 6)
    if (do_test) {
      // Find the vector perpendicular to the segment (p01 <-> p12)
      var perp : Vec3 = Vec3.cross(
        Util.trinorm(vertices[v0],vertices[v1],vertices[v2]),
        (p01 - p12)
      );
      // We want tmax1 to point opposite to perp, and
      // tmax0 and tmax2 to point along it.  Otherwise, exit out.
      if (Vec3.dot(tmax0 , perp) <= 0.0 ||
          Vec3.dot(tmax1 , perp) >= 0.0 ||
          Vec3.dot(tmax2 , perp) <= 0.0)
        return null;
    }
    // Fade lines
    k01 /= (k01 + thresh);
    k12 /= (k12 + thresh);
    return new Ridge(p01,k01,p12,k12);
  }

  // Draw apparent ridges of the mesh
  public function facesApparentRidges(
      ndotv : Array<Float>,
      t1q1 : Array<Vec3>, Dt1q1 : Array<Float>,
      do_bfcull : Bool, do_test : Bool, thresh : Float
  ) : Array<Ridge> {
    var ridges : Array<Ridge> =[];
    for (f in faces){
      // Draw apparent ridges in a triangle
      var v0 : Int = f[0];
      var v1 : Int = f[1];
      var v2 : Int = f[2];
      // Backface culling is turned off: getting contours from the
      // apparent ridge definition requires us to process faces that
      // may be (just barely) backfacing...
      if (do_bfcull &&
          ndotv[v0] <= 0 && ndotv[v1] <= 0 && ndotv[v2] <= 0){
        continue;
      }

      // Trivial reject if this face isn't getting past the threshold anyway
      var kmax0 : Float = t1q1[v0].z;
      var kmax1 : Float = t1q1[v1].z;
      var kmax2 : Float = t1q1[v2].z;
      if (kmax0 <= thresh && kmax1 <= thresh && kmax2 <= thresh)
        continue;
      
      // The "tmax" are the principal directions of view-dependent curvature,
      // flipped to point in the direction in which the curvature
      // is increasing.
      var emax0 : Float = Dt1q1[v0];
      var emax1 : Float = Dt1q1[v1];
      var emax2 : Float = Dt1q1[v2];
      var world_t1_0 : Vec3 = pdir1[v0]*t1q1[v0][0] + pdir2[v0]*t1q1[v0][1];
      var world_t1_1 : Vec3 = pdir1[v1]*t1q1[v1][0] + pdir2[v1]*t1q1[v1][1];
      var world_t1_2 : Vec3 = pdir1[v2]*t1q1[v2][0] + pdir2[v2]*t1q1[v2][1];

      var tmax0 : Vec3 = world_t1_0 * Dt1q1[v0];
      var tmax1 : Vec3 = world_t1_1 * Dt1q1[v1];
      var tmax2 : Vec3 = world_t1_2 * Dt1q1[v2];

      // We have a "zero crossing" if the tmaxes along an edge
      // point in opposite directions
      var z01 : Bool = (Vec3.dot(tmax0,tmax1) <= 0.0);
      var z12 : Bool = (Vec3.dot(tmax1,tmax2) <= 0.0);
      var z20 : Bool = (Vec3.dot(tmax2,tmax0) <= 0.0);

      if ((z01?1:0) + (z12?1:0) + (z20?1:0) < 2)
        continue;

      // Draw line segment
      if (!z01) {
        var r : Ridge = segmentApparentRidge(v1, v2, v0,
                  emax1, emax2, emax0,
                  kmax1, kmax2, kmax0,
                  tmax1, tmax2, tmax0,
                  thresh, false, do_test);
        if (r != null) ridges.push(r);
      } else if (!z12) {
        var r : Ridge = segmentApparentRidge(v2, v0, v1,
                  emax2, emax0, emax1,
                  kmax2, kmax0, kmax1,
                  tmax2, tmax0, tmax1,
                  thresh, false, do_test);
        if (r != null) ridges.push(r);
      } else if (!z20) {
        var r : Ridge = segmentApparentRidge(v0, v1, v2,
                  emax0, emax1, emax2,
                  kmax0, kmax1, kmax2,
                  tmax0, tmax1, tmax2,
                  thresh, false, do_test);
        if (r != null) ridges.push(r);
      } else {
        // All three edges have crossings -- connect all to center
        var r0 : Ridge = segmentApparentRidge(v1, v2, v0,
                  emax1, emax2, emax0,
                  kmax1, kmax2, kmax0,
                  tmax1, tmax2, tmax0,
                  thresh, true, do_test);
        var r1 : Ridge = segmentApparentRidge(v2, v0, v1,
                  emax2, emax0, emax1,
                  kmax2, kmax0, kmax1,
                  tmax2, tmax0, tmax1,
                  thresh, true, do_test);
        var r2 : Ridge = segmentApparentRidge(v0, v1, v2,
                  emax0, emax1, emax2,
                  kmax0, kmax1, kmax2,
                  tmax0, tmax1, tmax2,
                  thresh, true, do_test);
        if (r0 != null) ridges.push(r0);
        if (r1 != null) ridges.push(r1);
        if (r2 != null) ridges.push(r2);
      }
    }
    return ridges;
    
  }

  public function apparentRidges(
      eye : Vec3, thresh : Float
    ) : Array<Ridge> {
    var nv : Int = vertices.length;

    for (i in 0...nv){
      // Compute n DOT v
		  var viewdir : Vec3 = eye - vertices[i];
		  var rlv : Float = 1.0 / viewdir.len();
      viewdir *= rlv;
      ndotv[i] = Vec3.dot(viewdir,normals[i]);

      var u : Float = Vec3.dot(viewdir,pdir1[i]); 
      var u2 : Float = u*u;
      var v : Float = Vec3.dot(viewdir,pdir2[i]); 
      var v2 : Float = v*v;

      var csc2theta : Float = 1.0 / (u2 + v2);
      t1q1[i] = computeVertViewDepCurv(i, ndotv[i],
        u2*csc2theta, u*v*csc2theta, v2*csc2theta
      );
    }
    for (i in 0...nv){
      Dt1q1[i] = computeVertDt1q1(i,ndotv[i],t1q1);
    }
    return facesApparentRidges(
      ndotv,t1q1,Dt1q1,false,true,
      thresh/Util.sq(featureSize)
    );
  }

  public function computeBVHTrivial(){
    bvh = new BVHTree(this,faces.length);
    bvh.build();
  }
  public function computeBVH(){
    bvh = new BVHTree(this);
    bvh.build();
  }

  public function visible(eye: Vec3, p : Vec3, tolerance : Float = 2) : Bool{
    var epsilon = bsphere.r/Math.sqrt(faces.length)*tolerance;
    var r = new Ray();
    var x = p - eye;
    r.d = x.copy();
    r.d.normalize();
    r.o = eye;
    r.tmin = 0;
    r.tmax = x.len()-epsilon;
    var h = r.hitBVH(bvh);
    return h == null;
  }

}

@:expose
class Ray {
  public var o : Vec3;
  public var d : Vec3;
  public var tmin : Float;
  public var tmax : Float;
  public inline function new(){
  }
  public inline function hitBBox(bb:BBox):RayHit{
    var tx1 : Float = ((bb.min[0]) - (o[0])) / (d[0]);
    var tx2 : Float = ((bb.max[0]) - (o[0])) / (d[0]);
    var ty1 : Float = ((bb.min[1]) - (o[1])) / (d[1]);
    var ty2 : Float = ((bb.max[1]) - (o[1])) / (d[1]);
    var tz1 : Float = ((bb.min[2]) - (o[2])) / (d[2]);
    var tz2 : Float = ((bb.max[2]) - (o[2])) / (d[2]);
    
    var t1 = Math.max(Math.max(Math.min(tx1, tx2), Math.min(ty1, ty2)), Math.min(tz1, tz2));
    var t2 = Math.min(Math.min(Math.max(tx1, tx2), Math.max(ty1, ty2)), Math.max(tz1, tz2));
    if (t2 - t1 < 0){
      return null;
    }
    if (t1 > tmax || t2 < tmin){
      return null;
    }
    var h : RayHit = new RayHit(t1);
    h.t2 = t2;
    return h;
  }

  public inline function hitTriangle(p0:Vec3,p1:Vec3,p2:Vec3):RayHit{
    var e1 = p1 - p0;
    var e2 = p2 - p0;
    var s = o - p0;
    inline function det(a : Vec3, b : Vec3, c : Vec3) : Float{
      return Vec3.dot(Vec3.cross(a,b),c);
    }
    var _d : Vec3 = d * (-1.0);
    var denom : Float = det(e1,e2,_d);
    if (denom == 0){
      return null;
    }
    var uvt : Vec3 = new Vec3(det(s,e2,_d),det(e1,s,_d),det(e1,e2,s)) * (1/denom);
    var u : Float = uvt.x;
    var v : Float = uvt.y;
    var t : Float = uvt.z;
    if (u < 0 || v < 0 || (1-u-v) < 0 || t < tmin || t > tmax){
      return null;
    }
    var h : RayHit = new RayHit(t);
    h.u = u;
    h.v = v;
    return h;
  }
  public inline function hitBVH(bvh:BVHTree):RayHit{
    function hitNode(node:BVHNode):RayHit{
      if (node.isLeaf()){
        // return hitBBox(node.bbox);
        var tmin : Float = Math.POSITIVE_INFINITY;
        var closest : RayHit = null;
        for (i in node.begin...node.end){
          var h : RayHit = hitTriangle(
            bvh.mesh.vertices[bvh.faces[i][0]],
            bvh.mesh.vertices[bvh.faces[i][1]],
            bvh.mesh.vertices[bvh.faces[i][2]]
          );
          if (h != null){
            h.face = bvh.faces[i];
            if (tmin>h.t){
              tmin = h.t;
              closest = h;
            }
          }
        }
        return closest;
      }
      var hitL : RayHit = hitBBox(node.left.bbox);
      var hitR : RayHit = hitBBox(node.right.bbox);
      if (hitL != null && hitR == null){
        return hitNode(node.left);
      }else if (hitL == null && hitR != null){
        return hitNode(node.right);
      }else if (hitL == null && hitR == null){
        return null;
      }
      var first : BVHNode;
      var second : BVHNode;
      if (hitL.t < hitR.t){
        first = node.left;
        second = node.right;
      }else{
        first = node.right;
        second = node.left;
      }
      var h : RayHit = hitNode(first);
      if (h == null || h.t >= Math.max(hitL.t,hitR.t)){
        var h2 : RayHit = hitNode(second);
        if (h2 != null){
          if (h == null || h2.t < h.t){
            return h2;
          }
        }
      }
      return h;
    }
    return hitNode(bvh.root);
  }
}

@:expose
class RayHit {
  public var t : Float;
  public var t2 : Float;
  public var u : Float;
  public var v : Float;
  public var face : Face;
  public inline function new(_t:Float){
    t = _t;
  }
}

@:expose
class BBox {
  public var min : Vec3;
  public var max : Vec3;
  public inline function new(){
    min=new Vec3(Math.POSITIVE_INFINITY,Math.POSITIVE_INFINITY,Math.POSITIVE_INFINITY);
    max=new Vec3(Math.NEGATIVE_INFINITY,Math.NEGATIVE_INFINITY,Math.NEGATIVE_INFINITY);
  }
  public inline function centroid() : Vec3{
    return (min+max)*0.5;
  }
  public inline function add(p:Vec3){
    min.x = Math.min(min.x,p.x);
    min.y = Math.min(min.y,p.y);
    min.z = Math.min(min.z,p.z);
    max.x = Math.max(max.x,p.x);
    max.y = Math.max(max.y,p.y);
    max.z = Math.max(max.z,p.z);
  }
  public inline function merge(bb : BBox){
    min.x = Math.min(min.x,bb.min.x);
    min.y = Math.min(min.y,bb.min.y);
    min.z = Math.min(min.z,bb.min.z);
    max.x = Math.max(max.x,bb.max.x);
    max.y = Math.max(max.y,bb.max.y);
    max.z = Math.max(max.z,bb.max.z);
  }
  public inline function surfaceArea() : Float{
    var extent : Vec3 = max-min;
    var x : Float = Math.max(extent.x,0);
    var y : Float = Math.max(extent.y,0);
    var z : Float = Math.max(extent.z,0);
    return 2 * (x * z + x * y + y * z);
  }
}

@:expose
class BVHNode {
  public var left : BVHNode;
  public var right : BVHNode;
  public var begin : Int;
  public var end : Int;
  public var bbox : BBox;
  public inline function new (box:BBox,i0:Int,i1:Int){
    bbox = box;
    begin = i0;
    end = i1;
    left = null;
    right = null;
  }
  public inline function isLeaf(){
    return left == null && right == null;
  }
}

@:expose
class BVHTree {
  public var root : BVHNode;
  public var mesh : Mesh;
  public var faces : Array<Face>;
  public var maxLeafSize : Int;
  public var bucketCount : Int;

  public function new(
    _mesh : Mesh,
    _maxLeafSize:Int=4,_bucketCount:Int=8
  ){
    maxLeafSize = _maxLeafSize;
    bucketCount = _bucketCount;
    faces =_mesh.faces.slice(0); // shallow copy
    mesh = _mesh;
  }

  public function build(){
    function bboxAddFace(bbox:BBox,f : Face){
      for (j in 0...3){
        bbox.add(mesh.vertices[f[j]]);
      }
    }
    function buildRange(i0:Int,i1:Int) : BVHNode{
      var bbox:BBox = new BBox();
      for (i in i0...i1){
        bboxAddFace(bbox,faces[i]);
      }
      
      var node = new BVHNode(bbox,i0,i1);
      if (i1-i0<=maxLeafSize){
        return node;
      }
      var parts : Array<BVHPartition> = [];
      for (ax in 0...3){
        var buckets : Array<BVHBucket> = [];
        var lo : Float = bbox.min[ax];
        var hi : Float = bbox.max[ax];
        for (i in 0...bucketCount){
          var b : BVHBucket = new BVHBucket();
          b.min = lo + i / bucketCount * (hi - lo);
          b.max = (b.min) + (hi - lo) / bucketCount;
          buckets.push(b);
        }
        for (i in i0...i1){
          var bb = new BBox();
          bboxAddFace(bb,faces[i]);
          var c = bb.centroid();
          for (j in 0...bucketCount){
            if (buckets[j].min <= c[ax] && c[ax] <= buckets[j].max){
              buckets[j].count ++;
              buckets[j].bbox.merge(bb);
              buckets[j].area = buckets[j].area + bb.surfaceArea();
              break;
            }
          }
        }
        for (i in 0...bucketCount){
          var part : BVHPartition = new BVHPartition();
          part.planeIndex = i;
          part.axis = ax;
          for (j in 0...i){
            part.leftCount += buckets[j].count;
            part.leftArea += buckets[j].area;
            part.leftBBox.merge(buckets[j].bbox);
          }
          for (j in i...bucketCount){
            part.rightCount += buckets[j].count;
            part.rightArea += buckets[j].area;
            part.rightBBox.merge(buckets[j].bbox);
          }
          if (part.leftCount > 0 && part.rightCount > 0){
            part.SAH = part.leftBBox.surfaceArea()/part.leftCount 
                     + part.rightBBox.surfaceArea()/part.rightCount;
            parts.push(part);
          }
        }
      }
      if (parts.length == 0){
        return node;
      }
      var minSAH : Float = Math.POSITIVE_INFINITY;
      var minPart : BVHPartition = null;
      for (p in parts){
        if (p.SAH < minSAH){
          minSAH = p.SAH;
          minPart = p;
        }
      }
      function comp(f0:Face,f1:Face) : Int{
        var bb0 = new BBox();
        var bb1 = new BBox();
        bboxAddFace(bb0,f0);
        bboxAddFace(bb1,f1);
        var v = bb0.centroid()[minPart.axis]-bb1.centroid()[minPart.axis];
        if (v < 0) return -1;
        if (v > 0) return 1;
        return 0;
      }
      var sorted : Array<Face> = faces.slice(i0,i1);
      sorted.sort(comp);
      for (i in i0...i1){
        faces[i] = sorted[i-i0];
      }
      var m : Int = i0 + minPart.leftCount;
      node.left = buildRange(i0,m);
      node.right = buildRange(m,i1);
      
      return node;
    }
    root = buildRange(0,faces.length);
  }

}

class BVHBucket {
  public var min : Float;
  public var max : Float;
  public var count : Int;
  public var area : Float;
  public var bbox : BBox;
  public function new (){
    bbox = new BBox();
    area = 0;
    count = 0;
  }

}
class BVHPartition {
  public var planeIndex : Int;
  public var axis : Int;
  public var leftCount : Int = 0;
  public var rightCount : Int = 0;
  public var leftArea : Float = 0;
  public var rightArea : Float = 0;
  public var leftBBox : BBox;
  public var rightBBox : BBox;
  public var SAH : Float = 0;
  public function new(){
    leftBBox = new BBox();
    rightBBox = new BBox();
  }
}



// bare minimum .obj format parser
@:expose
class OBJParser {
#if sys
  public static function fromFile(path : String) : Mesh{
    return fromString(sys.io.File.getContent(path));
  }
#end
  public static function fromString(str : String) : Mesh{
    var mesh : Mesh = new Mesh();
    mesh.vertices = [];
    mesh.faces = [];

    var lines : Array<String> = str.split("\n");
    for (i in 0...lines.length){
      lines[i] = StringTools.trim(lines[i]);
      if (lines[i].charAt(0) == "#"){
        continue;
      }
      if (lines[i].length <= 2){
        continue;
      }
      var tok : Array<String> = lines[i].split(" ");
      var cmd : String = tok[0];
      if (cmd == "v"){
        var v : Vec3 = new Vec3(
          Std.parseFloat(tok[1]),
          Std.parseFloat(tok[2]),
          Std.parseFloat(tok[3])
        );
        mesh.vertices.push(v);
      }else if (cmd == "f"){
        var a : Int = Std.parseInt(tok[1].split("/")[0]);
        var b : Int = Std.parseInt(tok[2].split("/")[0]);
        var c : Int = Std.parseInt(tok[3].split("/")[0]);
        var nv : Int = mesh.vertices.length;
        mesh.faces.push(new Face(
          (a<0)?(nv+a):(a-1),
          (b<0)?(nv+b):(b-1),
          (c<0)?(nv+c):(c-1)
        ));
      }
    }
    return mesh;
  }
  
}




@:expose
class Line {
  public var x1 : Float;
  public var y1 : Float;
  public var x2 : Float;
  public var y2 : Float;
  public var opacity1 : Float = 1;
  public var opacity2 : Float = 1;
  public function new (_x1:Float, _y1:Float, _x2:Float, _y2:Float){
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
  }
  public function setOpacity(o1:Float,o2:Float){
    opacity1 = o1;
    opacity2 = o2;
  }
  public function flip(){
    var tmp : Float;
    tmp = x1;
    x1 = x2;
    x2 = tmp;
    tmp = y1;
    y1 = y2;
    y2 = tmp;
    tmp = opacity1;
    opacity1 = opacity2;
    opacity2 = tmp;
  }
}

@:expose
abstract Polyline (Array<Vec3>) {
  public var length(get, set):Int;
  public inline function new(){
    this = [];
  }
  public inline function get_length(){
    return this.length;
  }
  public inline function set_length(v:Int):Int{
    return this.length; //readonly
  }
  public inline function startY() : Int{
    return Util.rndInt(this[0].y);
  }
  public inline function endY() : Int{
    return Util.rndInt(this[this.length-1].y);
  }
  public inline function startX() : Float{
    return this[0].x;
  }
  public inline function endX() : Float{
    return this[this.length-1].x;
  }
  @:arrayAccess
  public inline function get(i:Int) {
    return this[i];
  }
  @:arrayAccess
  public inline function set(i:Int,v:Vec3):Vec3 {
    this[i] = v;
    return v;
  }
  public inline function push(v:Vec3):Int{
    this.push(v);
    return this.length;
  }
  public inline function unshift(v:Vec3){
    this.unshift(v);
  }
}


@:expose
class Render {
  public var mesh : Mesh;
  public var lines : Array<Line>; 
  public var polylines : Array<Polyline>;

  public var focal : Float = 1000;
  public var width : Int;
  public var height : Int;
  public var verbose : Bool = true;
 
  public var didPrecompute : Bool = false;

  public function new(_mesh:Mesh, w : Int, h : Int){
    mesh = _mesh;
    lines = [];
    width = w;
    height = h;
  }
  public function clear(){
    if (lines != null){
      lines.splice(0,lines.length);
    }
    if (polylines != null){
      polylines.splice(0,polylines.length);
    }
  }
  public function setFocal(f : Float){
    focal = f;
  }
  public function setVerbose(v : Int){
    verbose = (v > 0);
  }
  public function transform(mat4x4 : Array<Float>){
    for (i in 0...mesh.vertices.length){
      mesh.vertices[i] = Util.matTrfm(mat4x4,mesh.vertices[i]);
    }
  }
  public function scaleRotateTranslate(
    sx : Float, sy : Float, sz : Float,
    rx : Float, ry : Float, rz : Float,
    dx : Float, dy : Float, dz : Float
  ){
    var scl : Array<Float> = Util.matScal(sx,sy,sz);
    var rotx : Array<Float> = Util.matRotx(rx);
    var roty : Array<Float> = Util.matRoty(ry);
    var rotz : Array<Float> = Util.matRotz(rz);
    var trsl : Array<Float> = Util.matTrsl(dx,dy,dz);

    transform(scl);
    transform(
      Util.matMult(trsl,Util.matMult(rotz,Util.matMult(roty,rotx)))
    );
  }
  public function autoPlace(zFactor:Float=1.5,fFactor:Float=1.25){
    mesh.computeBSphere();
    transform(Util.matTrsl(
      -mesh.bsphere.o.x,-mesh.bsphere.o.y,-mesh.bsphere.o.z
    ));
    var r = Util.min(width,height)/2;
    // trace(mesh.bsphere.o,mesh.bsphere.r);
    transform(Util.matScal(
      r/mesh.bsphere.r,r/mesh.bsphere.r,r/mesh.bsphere.r
    ));
    transform(Util.matTrsl(0,0,r*zFactor));

    // mesh.computeBSphere();
    // trace(mesh.bsphere.o,mesh.bsphere.r);
    setFocal(r*fFactor);
  }
  public function vertices(){
    var offs = new Vec3(width/2,height/2,0);
    var yflip = new Vec3(-1,-1,1);
    for (i in 0...mesh.vertices.length){
      var v : Vec3 = mesh.vertices[i];
      var p : Vec3 = Util.matProj(focal,v)*yflip+offs;
      lines.push(new Line(p.x-1,p.y-1,p.x+1,p.y+1));
      lines.push(new Line(p.x+1,p.y-1,p.x-1,p.y+1));
    }
  }
  public function edges(){
    var offs = new Vec3(width/2,height/2,0);
    var yflip = new Vec3(-1,-1,1);
    for (i in 0...mesh.faces.length){
      var f : Face = mesh.faces[i];

      var p0 : Vec3 = Util.matProj(focal,mesh.vertices[f[0]])*yflip+offs;
      var p1 : Vec3 = Util.matProj(focal,mesh.vertices[f[1]])*yflip+offs;
      var p2 : Vec3 = Util.matProj(focal,mesh.vertices[f[2]])*yflip+offs;
      
      lines.push(new Line(p0.x,p0.y,p1.x,p1.y));
      lines.push(new Line(p1.x,p1.y,p2.x,p2.y));
      lines.push(new Line(p2.x,p2.y,p0.x,p0.y));
    }
  }
  public function apparentRidges(thresh : Float, cull : Float = 2){
    if (!didPrecompute){
      if (verbose)trace("precomputing mesh properties...");
      mesh.precompute(cull>=0,verbose);
      didPrecompute = true;
    }
    var offs = new Vec3(width/2,height/2,0);
    var yflip = new Vec3(-1,-1,1);
    var eye = new Vec3(0,0,0);

    if (verbose)trace("generating apparent ridges...");
    var ridges : Array<Ridge> = mesh.apparentRidges(eye,thresh);
    if (verbose)trace("projecting apparent ridges onto 2D plane...");
    for (i in 0...ridges.length){
      if (cull >= 0){
        if (!mesh.visible(eye,ridges[i].A,cull) && !mesh.visible(eye,ridges[i].B,cull)){
          continue;
        }
      }
      var p0 : Vec3 = Util.matProj(focal,ridges[i].A)*yflip+offs;
      var p1 : Vec3 = Util.matProj(focal,ridges[i].B)*yflip+offs;
      var l : Line = new Line(p0.x,p0.y,p1.x,p1.y);
      l.opacity1 = ridges[i].strengthA;
      l.opacity2 = ridges[i].strengthB;
      lines.push(l);
    }
    if (verbose)trace("apparent ridges computation finished.");
  }

  public function buildPolylines(epsilon : Float = 1){
    if (verbose)trace("building polylines from ridge segments...");
    polylines = [];
    lines = lines.filter(function(a:Line):Bool{
      return !(
        a.y1 < 0 || a.y1 > height-1 || 
        a.y2 < 0 || a.y2 > height-1
      );
    });

    for (i in 0...lines.length){
      var y1 : Float = lines[i].y1;
      var y2 : Float = lines[i].y2;
      if (y1 > y2){
        lines[i].flip();
      }else if (y1 == y2){
        if (lines[i].x1 > lines[i].x2){
          lines[i].flip();
        }
      }
    }
    var rows : Array<Array<Polyline>> = [for (i in 0...height) []];
    var ends : Array<Array<Polyline>> = [for (i in 0...height) []];

    function singleton(a:Line):Polyline{
      var p : Polyline = new Polyline();
      p.push(new Vec3(a.x1,a.y1,a.opacity1));
      p.push(new Vec3(a.x2,a.y2,a.opacity2));
      return p;
    }

    for (i in 0...lines.length){
      var p : Polyline = singleton(lines[i]);
      rows[Util.rndInt(lines[i].y1)].push(p);
      ends[Util.rndInt(lines[i].y2)].push(p);
    }
    for (i in 0...rows.length){
      var nj : Int = rows[i].length;
      for (_j in 0...nj){    var j : Int = nj-_j-1;
        if (rows[i][j] == null){
          continue;
        }
        var nk : Int = ends[i].length;
        for (_k in 0...nk){  var k : Int = nk-_k-1;
          if (ends[i][k] == null){
            continue;
          }

          if (rows[i][j] == ends[i][k]){
            continue;
          }
          var r : Int = rows[i][j].endY();
          var d : Float = Math.abs( rows[i][j].startX() - ends[i][k].endX() );
          if ( d <= epsilon){
            if (d < 1){
              ends[i][k][ends[i][k].length-1].z=(
                ends[i][k][ends[i][k].length-1].z+
                rows[i][j][0].z
              )/2;
            }
            for (t in (d<1?1:0)...rows[i][j].length){
              ends[i][k].push(rows[i][j][t]);
            }
            ends[r].remove(rows[i][j]);
            ends[r].push(ends[i][k]);
            ends[i][k] = null;
            
            break;
          }
        }
      }
    }
    for (i in 0...ends.length){
      for (j in 0...ends[i].length){
        if (ends[i][j] != null){
          polylines.push(ends[i][j]);
        }
      }
    }
    polylines = polylines.filter(function(p:Polyline):Bool{
      if (p.length > 2){
        return true;
      }
      if (p.length < 2){
        return false;
      }
      if (Vec3.dist(p[0],p[1])<epsilon){
        return false;
      }
      return true;
    });

    if (verbose)trace("polylines built.");
  }



}

@:expose
class PixelMap {

  public static function raycast(render:Render,fun:RayHit->Int->Int->Void){
    var min : Float = Math.POSITIVE_INFINITY;
    var max : Float = Math.NEGATIVE_INFINITY;
    var width  : Int = render.width;
    var height : Int  = render.height;
    var hw : Int = Std.int(width/2);
    var hh : Int = Std.int(height/2);

    for (y in -hh...(height-hh)){
      for (x in -hw...(width-hw)){
        var r : Ray = new Ray();
        r.o = new Vec3(0,0,0);
        r.d = new Vec3(-x,-y,render.focal);
        r.tmax = Math.POSITIVE_INFINITY;
        r.tmin = 0;
        r.d.normalize();

        var h : RayHit = r.hitBVH(render.mesh.bvh);

        fun(h,x+hh,y+hh);
      }
    }
  }
  
  public static function depth(render:Render,normalize:Bool=false) : haxe.ds.Vector<Float>{
    var data : haxe.ds.Vector<Float> = new haxe.ds.Vector<Float>(render.width*render.height);

    var min : Float = Math.POSITIVE_INFINITY;
    var max : Float = Math.NEGATIVE_INFINITY;

    
    raycast(render,function (h:RayHit,x:Int,y:Int){
      if (h == null){
        data[y*render.width+x]=Math.POSITIVE_INFINITY;
      }else{
        min = Math.min(min,h.t);
        max = Math.max(max,h.t);
        data[y*render.width+x] = h.t;
      }
    });

    if (normalize){
      for (i in 0...data.length){
        if (data[i] != Math.POSITIVE_INFINITY){
          data[i] = 1-(data[i]-min)/(max-min);
        }else{
          data[i] = 0;
        }
      }
    }
    return data;
  }

  public static function normal(render:Render) : haxe.ds.Vector<Float>{
    var data : haxe.ds.Vector<Float> = new haxe.ds.Vector<Float>(render.width*render.height*3);

    raycast(render,function (h:RayHit,x:Int,y:Int){
      var idx : Int = (y*render.width+x)*3;
      if (h == null){
        data[idx  ]=0;
        data[idx+1]=0;
        data[idx+2]=0;
      }else{
        var f : Face = h.face;
        var n0 : Vec3 = render.mesh.normals[f[0]];
        var n1 : Vec3 = render.mesh.normals[f[1]];
        var n2 : Vec3 = render.mesh.normals[f[2]];
        var n : Vec3 = n0 * (1-h.u-h.v) + n1 * h.u + n2 * h.v;
        data[idx  ] = n[0];
        data[idx+1] = n[1];
        data[idx+2] = n[2];
      }
    });

    return data;
  }

  
  public static function curvature(render:Render) : haxe.ds.Vector<Float>{
    var data : haxe.ds.Vector<Float> = new haxe.ds.Vector<Float>(render.width*render.height*2);

    raycast(render,function (h:RayHit,x:Int,y:Int){
      var idx : Int = (y*render.width+x)*2;
      if (h == null){
        data[idx  ]=0;
        data[idx+1]=0;
      }else{
        var f : Face = h.face;
        var c1a : Float = render.mesh.curv1[f[0]];
        var c1b : Float = render.mesh.curv1[f[1]];
        var c1c : Float = render.mesh.curv1[f[2]];
        var c2a : Float = render.mesh.curv2[f[0]];
        var c2b : Float = render.mesh.curv2[f[1]];
        var c2c : Float = render.mesh.curv2[f[2]];
        var c1 : Float = c1a * (1-h.u-h.v) + c1b * h.u + c1c * h.v;
        var c2 : Float = c2a * (1-h.u-h.v) + c2b * h.u + c2c * h.v;
        data[idx  ] = c1;
        data[idx+1] = c2;
      }
    });

    return data;
  }

  public static function lambertian(render:Render,light:Vec3,normalize:Bool=true) : haxe.ds.Vector<Float>{
    var data : haxe.ds.Vector<Float> = new haxe.ds.Vector<Float>(render.width*render.height);
    var min : Float = Math.POSITIVE_INFINITY;
    var max : Float = Math.NEGATIVE_INFINITY;

    raycast(render,function (h:RayHit,x:Int,y:Int){
      var idx : Int = (y*render.width+x);
      if (h == null){
        data[idx]=Math.NEGATIVE_INFINITY;
      }else{
        var f : Face = h.face;
        var n0 : Vec3 = render.mesh.normals[f[0]];
        var n1 : Vec3 = render.mesh.normals[f[1]];
        var n2 : Vec3 = render.mesh.normals[f[2]];
        var n : Vec3 = n0 * (1-h.u-h.v) + n1 * h.u + n2 * h.v;
        var ndotl : Float = Vec3.dot(n,light);
        min = Math.min(min,ndotl);
        max = Math.max(max,ndotl);
        data[idx] = ndotl;
      }
    });
    if (normalize){
      for (i in 0...data.length){
        if (data[i] != Math.NEGATIVE_INFINITY){
          data[i] = (data[i]-min)/(max-min);
        }else{
          data[i] = 0;
        }
      }
    }
    return data;
  }

  public static function ambientOcclusion(render:Render,numSamples:Int=32,normalize:Bool=true) : haxe.ds.Vector<Float>{
    var data : haxe.ds.Vector<Float> = new haxe.ds.Vector<Float>(render.width*render.height);
    var min : Float = Math.POSITIVE_INFINITY;
    var max : Float = Math.NEGATIVE_INFINITY;

    raycast(render,function (h:RayHit,x:Int,y:Int){
      var idx : Int = (y*render.width+x);
      if (h == null){
        data[idx]=Math.NEGATIVE_INFINITY;
      }else{
        var f : Face = h.face;

        var p0 : Vec3 = render.mesh.vertices[f[0]];
        var p1 : Vec3 = render.mesh.vertices[f[1]];
        var p2 : Vec3 = render.mesh.vertices[f[2]];
        var o : Vec3 = p0 * (1-h.u-h.v) + p1 * h.u + p2 * h.v;
        
        var cnt : Int = 0;
        for (i in 0...numSamples){
          var d : Vec3 = Util.uniformHemisphereSampler();
          if (Math.random()<0.5){
            d.z = -d.z;
          }
          var r : Ray = new Ray();
          r.d = d;
          r.o = o;
          r.tmin = 0.1;
          r.tmax = Math.POSITIVE_INFINITY;
          var h = r.hitBVH(render.mesh.bvh);
          if (h != null){
            cnt ++;
          }
        }
        var v : Float = 1-cnt/numSamples;
        min = Math.min(min,v);
        max = Math.max(max,v);
        data[idx] = v;
      }
    });

    if (normalize){
      for (i in 0...data.length){
        if (data[i] != Math.NEGATIVE_INFINITY){
          data[i] = (data[i]-min)/(max-min);
        }else{
          data[i] = 0;
        }
      }
    }
    return data;
  }


  public static function toPPMString(
    data : haxe.ds.Vector<Float>, 
    w : Int, h : Int, min : Float, max : Float
  ){
    inline function nor2int(x : Float){
      return Std.int( Math.min(Math.max(((x-min)/(max-min))*255,0),255) );
    }
    var chan : Int = Std.int(data.length/(w*h));

    var out : StringBuf = new StringBuf();
    out.add('P3\n$w $h\n255\n');
    for (i in 0...Std.int(data.length/chan)){
      var u : Int = nor2int(data[i*chan]);
      var v : Int = (chan==1)?u:nor2int(data[i*chan+1]);
      var w : Int = (chan==1)?u:((chan==2)?128:nor2int(data[i*chan+2]));
      out.add('$u $v $w ');
    }
    return out.toString();
  }
}

@:expose
class SVGWriter {

  static inline function rd(x : Float){
    return Math.round(x*100)/100;
  }

  public static function lines(render : Render, useOpacity : Bool = true) : String{
    var w : Int = render.width;
    var h : Int = render.height;
    var out : StringBuf = new StringBuf();
    out.add('<svg version="1.1" xmlns="http://www.w3.org/2000/svg" width="$w" height="$h">\n');

    for (i in 0...render.lines.length){
      var x1 : Float = rd(render.lines[i].x1);
      var y1 : Float = rd(render.lines[i].y1);
      var x2 : Float = rd(render.lines[i].x2);
      var y2 : Float = rd(render.lines[i].y2);
      var o : Float = useOpacity ? ((render.lines[i].opacity1+render.lines[i].opacity2)/2) : 1;
      var oi : Int = Std.int(255-o*255);
      out.add('  <line x1="$x1" y1="$y1" x2="$x2" y2="$y2" fill="none" stroke="rgb($oi,$oi,$oi)" stroke-width="1" stroke-linecap="round"/>\n');
    }
    out.add("</svg>\n");
    return out.toString();
  }

  public static function polylines(render : Render, colorful : Bool = false) : String{
    var w : Int = render.width;
    var h : Int = render.height;
    var out : StringBuf = new StringBuf();
    out.add('<svg version="1.1" xmlns="http://www.w3.org/2000/svg" width="$w" height="$h">\n');

    var color : String = "black";
    for (i in 0...render.polylines.length){
      out.add('  <polyline points="');
      for (j in 0...render.polylines[i].length){
        var p : Vec3 = render.polylines[i][j];
        out.add(""+rd(p.x)+","+rd(p.y)+" ");
      }
      if (colorful){
        color = 'rgb(${Std.int(Math.random()*128)},${Std.int(Math.random()*128)},${Std.int(Math.random()*128)})';
      }
      out.add('" fill="none" stroke="$color" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>\n');
    }
    out.add("</svg>\n");
    return out.toString();
  }

  public static function gradients(render : Render, acc : Float = 1) : String{
    var w : Int = render.width;
    var h : Int = render.height;
    var out : StringBuf = new StringBuf();
    out.add('<svg version="1.1" xmlns="http://www.w3.org/2000/svg" width="$w" height="$h">\n');
    var gid : Int = 0;
    for (i in 0...render.polylines.length){
      for (j in 0...render.polylines[i].length){
        var p : Vec3 = render.polylines[i][j];
        var oi : Int = Std.int(255-p.z*255);
        out.add('  <circle cx="${rd(p.x)}" cy="${rd(p.y)}" r="0.5" stroke="none" fill="rgb($oi,$oi,$oi)"/>\n');
      }
    }
    for (i in 0...render.polylines.length){
      out.add('  <g fill="none" stroke-width="1">\n');
      for (j in 0...render.polylines[i].length-1){
        var p : Vec3 = render.polylines[i][j];
        var q : Vec3 = render.polylines[i][j+1];
        var d : Int = Std.int(Math.abs(p.z - q.z)*10*acc)+1;

        for (k in 0...d){
          var t : Float = k/d-0.01;
          var s : Float = (k+1)/d+0.01;
          var a : Vec3 = p * (1-t) + q * t;
          var b : Vec3 = p * (1-s) + q * s;
          var o : Float = (a.z + b.z)/2;
          var oi : Int = Std.int(255-o*255);
          out.add('    <line x1="${rd(a.x)}" y1="${rd(a.y)}" x2="${rd(b.x)}" y2="${rd(b.y)}" stroke="rgb($oi,$oi,$oi)"/>\n');
        }
        gid ++;
      }
      out.add('  </g>\n');
    }
    out.add("</svg>\n");
    return out.toString();
  }


}
