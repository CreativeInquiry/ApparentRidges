﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_apparentridges__OBJParser
#define INC_apparentridges__OBJParser
typedef struct _apparentridges__$OBJParser *apparentridges__$OBJParser;
typedef struct _apparentridges__OBJParser *apparentridges__OBJParser;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>
#include <apparentridges/Mesh.h>


struct _apparentridges__$OBJParser {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
	vclosure* fromFile;
	vclosure* fromString;
};
struct _apparentridges__OBJParser {
	hl_type *$type;
};
#endif

