﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_apparentridges__RayHit
#define INC_apparentridges__RayHit
typedef struct _apparentridges__$RayHit *apparentridges__$RayHit;
typedef struct _apparentridges__RayHit *apparentridges__RayHit;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>
#include <hl/types/ArrayBytes_Int.h>


struct _apparentridges__$RayHit {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
};
struct _apparentridges__RayHit {
	hl_type *$type;
	double _t;
	double t2;
	double u;
	double v;
	hl__types__ArrayBytes_Int face;
};
#endif

