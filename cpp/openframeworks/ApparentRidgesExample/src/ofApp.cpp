
#include <hxcpp.h>
#include "hxMath.h"

#ifndef INCLUDED_Std
#include <Std.h>
#endif
#ifndef INCLUDED_StringBuf
#include <StringBuf.h>
#endif
#ifndef INCLUDED_apparentridges_Line
#include <apparentridges/Line.h>
#endif
#ifndef INCLUDED_apparentridges_Render
#include <apparentridges/Render.h>
#endif
#ifndef INCLUDED_apparentridges_SVGWriter
#include <apparentridges/SVGWriter.h>
#endif
#ifndef INCLUDED_apparentridges_Mesh
#include <apparentridges/Mesh.h>
#endif
#ifndef INCLUDED_apparentridges_OBJParser
#include <apparentridges/OBJParser.h>
#endif
#ifndef INCLUDED_apparentridges_Util
#include <apparentridges/Util.h>
#endif

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif

#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
  
  // HX_TOP_OF_STACK
  // ::hx::Boot();
	// __boot_all();
  hx::StackContext *_hx_ctx = hx::gMainThreadContext;

  // mesh = ::apparentridges::OBJParser_obj::fromFile(HX_("../../../data/lucy_100K.obj",ff,ff,ff,ff));
  // render = ::apparentridges::Render_obj::__alloc(HX_CTX,mesh,800,800);
  
  // render->focal = 400;
  // render->transform(::apparentridges::Util_obj::matTrsl(0,-80,0));
  // render->transform(::apparentridges::Util_obj::matScal(4,4,4));
  // // render->transform(::apparentridges::Util_obj::matRoty(1.57));
  // render->setVerbose(0);
  

  mesh = ::apparentridges::OBJParser_obj::fromFile(HX_("../../../data/nefertiti_100K.obj",ff,ff,ff,ff));
  render = ::apparentridges::Render_obj::__alloc(HX_CTX,mesh,800,800);
  
  render->focal = 400;
  render->transform(::apparentridges::Util_obj::matRotx(-M_PI*0.5));
  render->transform(::apparentridges::Util_obj::matRoty(M_PI*1.1));
  // render->transform(::apparentridges::Util_obj::matTrsl(0,-80,0));
  // render->transform(::apparentridges::Util_obj::matScal(4,4,4));
  // render->transform(::apparentridges::Util_obj::matRoty(1.57));
  render->setVerbose(0);
  

	gui.setup();
	gui.add(thresh.setup("threshold", 0.5, 0, 1));

}

//--------------------------------------------------------------
void ofApp::update(){

  
}

//--------------------------------------------------------------
void ofApp::draw(){
  
  ofBackground(255);
  ofNoFill();
  ofSetLineWidth(1);

  ofSetColor(0,0,0);

  render->transform(::apparentridges::Util_obj::matRoty(0.1));
  render->transform(::apparentridges::Util_obj::matTrsl(0,0,500));

  // render->didPrecompute = false;
  if (render->didPrecompute){
    mesh->computeNormals();
    mesh->computeCurvatures();
    mesh->computeBVHTrivial();
  }
  render->apparentRidges(thresh,-1);
  double x1 = render->lines->__get(0).StaticCast<  ::apparentridges::Line >()->x1;
  int len = render->lines->length;
  for (int i = 0; i < len; i++){
    double x1 = render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->x1 ;
    double y1 = render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->y1 ;
    double x2 = render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->x2 ;
    double y2 = render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->y2 ;
    // printf("%f %f %f %f\n",x1,y1,x2,y2);
    ofDrawLine(x1,y1,x2,y2);
  }
  render->transform(::apparentridges::Util_obj::matTrsl(0,0,-500));
  render->clear();

  gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
