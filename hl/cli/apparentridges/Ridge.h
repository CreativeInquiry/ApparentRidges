﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_apparentridges__Ridge
#define INC_apparentridges__Ridge
typedef struct _apparentridges__$Ridge *apparentridges__$Ridge;
typedef struct _apparentridges__Ridge *apparentridges__Ridge;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>
#include <hl/types/ArrayBytes_Float.h>


struct _apparentridges__$Ridge {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
};
struct _apparentridges__Ridge {
	hl_type *$type;
	hl__types__ArrayBytes_Float A;
	hl__types__ArrayBytes_Float B;
	double strengthA;
	double strengthB;
};
#endif

