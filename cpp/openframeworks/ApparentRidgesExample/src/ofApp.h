
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

#pragma once

#include "ofMain.h"
#include "ofxGui.h"

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

		::apparentridges::Mesh mesh;
		::apparentridges::Render render;

		ofxPanel gui;
		ofxFloatSlider thresh;
		
};
