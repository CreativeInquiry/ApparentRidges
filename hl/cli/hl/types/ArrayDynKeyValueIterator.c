﻿// Generated by HLC 4.1.3 (HL v4)
#define HLC_BOOT
#include <hlc.h>
#include <hl/types/ArrayDynKeyValueIterator.h>
extern hl_type t$hl_types_ArrayDyn;
extern hl_type t$_dyn;
void haxe_iterators_ArrayKeyValueIterator_new(haxe__iterators__ArrayKeyValueIterator,hl__types__ArrayDyn);

void hl_types_ArrayDynKeyValueIterator_new(hl__types__ArrayDynKeyValueIterator r0,hl__types__ArrayBase r1) {
	hl__types__ArrayDyn r4;
	vdynamic *r3;
	r3 = NULL;
	r4 = (hl__types__ArrayDyn)hl_dyn_castp(&r3,&t$_dyn,&t$hl_types_ArrayDyn);
	haxe_iterators_ArrayKeyValueIterator_new(((haxe__iterators__ArrayKeyValueIterator)r0),r4);
	r0->a = r1;
	return;
}

