﻿// Generated by HLC 4.1.3 (HL v4)
#ifndef INC_hl__types__ArrayDynKeyValueIterator
#define INC_hl__types__ArrayDynKeyValueIterator
typedef struct _hl__types__$ArrayDynKeyValueIterator *hl__types__$ArrayDynKeyValueIterator;
typedef struct _hl__types__ArrayDynKeyValueIterator *hl__types__ArrayDynKeyValueIterator;
#include <hl/Class.h>
#include <hl/BaseType.h>
#include <_std/String.h>
#include <haxe/iterators/ArrayKeyValueIterator.h>
#include <hl/types/ArrayDyn.h>
#include <hl/types/ArrayBase.h>


struct _hl__types__$ArrayDynKeyValueIterator {
	hl_type *$type;
	hl_type* __type__;
	vdynamic* __meta__;
	varray* __implementedBy__;
	String __name__;
	vdynamic* __constructor__;
};
struct _hl__types__ArrayDynKeyValueIterator {
	hl_type *$type;
	int current;
	hl__types__ArrayDyn array;
	hl__types__ArrayBase a;
};
#endif

