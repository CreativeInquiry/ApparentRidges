
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

#include "ofMain.h"
#include "ofApp.h"


//========================================================================
int main( ){
  HX_TOP_OF_STACK
  ::hx::Boot();
	__boot_all();

	ofSetupOpenGL(1024,768,OF_WINDOW);			// <-------- setup the GL context

	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp(new ofApp());

}
