﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_apparentridges__SVGWriter
#define INC_apparentridges__SVGWriter
typedef struct _apparentridges__$SVGWriter *apparentridges__$SVGWriter;
typedef struct _apparentridges__SVGWriter *apparentridges__SVGWriter;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>
#include <apparentridges/Render.h>


struct _apparentridges__$SVGWriter {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
	vclosure* rd;
	vclosure* lines;
	vclosure* polylines;
	vclosure* gradients;
};
struct _apparentridges__SVGWriter {
	hl_type *$type;
};
#endif

