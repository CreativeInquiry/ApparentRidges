﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_apparentridges__Line
#define INC_apparentridges__Line
typedef struct _apparentridges__$Line *apparentridges__$Line;
typedef struct _apparentridges__Line *apparentridges__Line;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>


struct _apparentridges__$Line {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
};
struct _apparentridges__Line {
	hl_type *$type;
	double x1;
	double y1;
	double x2;
	double y2;
	double opacity1;
	double opacity2;
};
#endif

