﻿// Generated by HLC 4.1.3 (HL v4)
#define HLC_BOOT
#include <hlc.h>
#include <_std/StringTools.h>
#include <hl/natives.h>
vdynamic* String_charCodeAt(String,int);
extern hl_type t$_i32;
String String_substr(String,int,vdynamic*);
#include <_std/StringBuf.h>
extern hl_type t$StringBuf;
void StringBuf_new(StringBuf);
void StringBuf_add(StringBuf,vdynamic*);
String StringBuf_toString(StringBuf);
#include <hl/types/ArrayObj.h>
hl__types__ArrayObj String_split(String,String);
String hl_types_ArrayObj_join(hl__types__ArrayObj,String);

bool StringTools_startsWith(String r0,String r1) {
	bool r4;
	vbyte *r5, *r6;
	int r2, r3, r7, r8;
	if( r0 == NULL ) hl_null_access();
	r2 = r0->length;
	if( r1 == NULL ) hl_null_access();
	r3 = r1->length;
	if( r2 < r3 ) goto label$de37252_1_19;
	r5 = r0->bytes;
	r2 = 0;
	r6 = r1->bytes;
	r3 = 0;
	r7 = r1->length;
	r8 = 1;
	r7 = r7 << r8;
	r2 = hl_bytes_compare(r5,r2,r6,r3,r7);
	r3 = 0;
	if( r2 == r3 ) goto label$de37252_1_17;
	r4 = false;
	goto label$de37252_1_18;
	label$de37252_1_17:
	r4 = true;
	label$de37252_1_18:
	return r4;
	label$de37252_1_19:
	r4 = false;
	return r4;
}

bool StringTools_isSpace(String r0,int r1) {
	bool r5;
	vdynamic *r2, *r6;
	int r3, r4;
	if( r0 == NULL ) hl_null_access();
	r2 = String_charCodeAt(r0,r1);
	r3 = r2 ? r2->v.i : 0;
	r4 = 8;
	if( r4 >= r3 ) goto label$de37252_2_8;
	r3 = r2 ? r2->v.i : 0;
	r4 = 14;
	if( r3 < r4 ) goto label$de37252_2_15;
	label$de37252_2_8:
	r3 = 32;
	r6 = hl_alloc_dynamic(&t$_i32);
	r6->v.i = r3;
	if( r2 == r6 || (r2 && r6 && (r2->v.i == r6->v.i)) ) goto label$de37252_2_13;
	r5 = false;
	goto label$de37252_2_14;
	label$de37252_2_13:
	r5 = true;
	label$de37252_2_14:
	return r5;
	label$de37252_2_15:
	r5 = true;
	return r5;
}

String StringTools_ltrim(String r0) {
	String r2;
	bool r5;
	vdynamic *r6;
	int r1, r3, r4;
	if( r0 == NULL ) hl_null_access();
	r1 = r0->length;
	r3 = 0;
	label$de37252_3_3:
	if( r3 >= r1 ) goto label$de37252_3_9;
	r5 = StringTools_isSpace(r0,r3);
	if( !r5 ) goto label$de37252_3_9;
	++r3;
	goto label$de37252_3_3;
	label$de37252_3_9:
	r4 = 0;
	if( r4 >= r3 ) goto label$de37252_3_16;
	if( r0 == NULL ) hl_null_access();
	r4 = r1 - r3;
	r6 = hl_alloc_dynamic(&t$_i32);
	r6->v.i = r4;
	r2 = String_substr(r0,r3,r6);
	return r2;
	label$de37252_3_16:
	return r0;
}

String StringTools_rtrim(String r0) {
	String r2;
	bool r6;
	vdynamic *r7;
	int r1, r3, r4, r5;
	if( r0 == NULL ) hl_null_access();
	r1 = r0->length;
	r3 = 0;
	label$de37252_4_3:
	if( r3 >= r1 ) goto label$de37252_4_12;
	r4 = r1 - r3;
	r5 = 1;
	r4 = r4 - r5;
	r6 = StringTools_isSpace(r0,r4);
	if( !r6 ) goto label$de37252_4_12;
	++r3;
	goto label$de37252_4_3;
	label$de37252_4_12:
	r5 = 0;
	if( r5 >= r3 ) goto label$de37252_4_20;
	if( r0 == NULL ) hl_null_access();
	r4 = 0;
	r5 = r1 - r3;
	r7 = hl_alloc_dynamic(&t$_i32);
	r7->v.i = r5;
	r2 = String_substr(r0,r4,r7);
	return r2;
	label$de37252_4_20:
	return r0;
}

String StringTools_trim(String r0) {
	String r1;
	r1 = StringTools_rtrim(r0);
	r1 = StringTools_ltrim(r1);
	return r1;
}

String StringTools_rpad(String r0,String r1,int r2) {
	String r5;
	StringBuf r7;
	int r4, r6;
	if( r1 == NULL ) hl_null_access();
	r4 = r1->length;
	r6 = 0;
	if( r6 < r4 ) goto label$de37252_6_5;
	return r0;
	label$de37252_6_5:
	r7 = (StringBuf)hl_alloc_obj(&t$StringBuf);
	StringBuf_new(r7);
	StringBuf_add(r7,((vdynamic*)r0));
	label$de37252_6_8:
	if( r7 == NULL ) hl_null_access();
	r4 = r7->pos;
	r6 = 1;
	r4 = r4 >> r6;
	if( r4 >= r2 ) goto label$de37252_6_16;
	StringBuf_add(r7,((vdynamic*)r1));
	goto label$de37252_6_8;
	label$de37252_6_16:
	r5 = StringBuf_toString(r7);
	return r5;
}

String StringTools_replace(String r0,String r1,String r2) {
	String r3;
	hl__types__ArrayObj r4;
	if( r0 == NULL ) hl_null_access();
	r4 = String_split(r0,r1);
	if( r4 == NULL ) hl_null_access();
	r3 = hl_types_ArrayObj_join(r4,r2);
	return r3;
}

