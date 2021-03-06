﻿// Generated by HLC 4.1.3 (HL v4)
#define HLC_BOOT
#include <hlc.h>
#include <hl/BaseType.h>
#include <hl/Class.h>
#include <_std/String.h>
#include <_std/Date.h>
#include <_std/Std.h>
#include <_std/Math.h>
#include <hl/types/ArrayAccess.h>
#include <hl/types/ArrayBase.h>
#include <hl/types/ArrayBytes_hl_F32.h>
#include <_std/StringBuf.h>
#include <_std/SysError.h>
#include <_std/Sys.h>
#include <hl/natives.h>
#include <hl/Enum.h>
#include <apparentridges/Ridge.h>
#include <apparentridges/BSphere.h>
#include <apparentridges/BBox.h>
#include <apparentridges/BVHNode.h>
#include <apparentridges/BVHTree.h>
#include <apparentridges/Mesh.h>
#include <haxe/Log.h>
#include <apparentridges/RayHit.h>
#include <apparentridges/Ray.h>
#include <apparentridges/BVHPartition.h>
#include <apparentridges/BVHBucket.h>
#include <apparentridges/Line.h>
#include <apparentridges/Render.h>
#include <haxe/Exception.h>
#include <haxe/ValueException.h>
#include <haxe/io/Encoding.h>
#include <haxe/io/Bytes.h>
#include <haxe/io/BytesDataImpl.h>
#include <haxe/io/Eof.h>
#include <haxe/io/Output.h>
#include <haxe/io/Error.h>
#include <haxe/iterators/ArrayIterator.h>
#include <haxe/iterators/ArrayKeyValueIterator.h>
#include <hl/NativeArrayIterator_Dynamic.h>
#include <hl/NativeArrayIterator_Int.h>
#include <hl/types/BytesIterator_Float.h>
#include <hl/types/BytesKeyValueIterator_Float.h>
#include <hl/types/BytesIterator_Int.h>
#include <hl/types/BytesKeyValueIterator_Int.h>
#include <hl/types/BytesIterator_hl_F32.h>
#include <hl/types/BytesKeyValueIterator_hl_F32.h>
#include <hl/types/BytesIterator_hl_UI16.h>
#include <hl/types/BytesKeyValueIterator_hl_UI16.h>
#include <hl/types/ArrayDynIterator.h>
#include <hl/types/ArrayDynKeyValueIterator.h>
#include <hl/types/ArrayObjIterator.h>
#include <hl/types/ArrayObjKeyValueIterator.h>
#include <sys/io/FileOutput.h>
#include <hl/CoreType.h>
#include <hl/CoreEnum.h>
#include <_std/StringTools.h>
#include <hl/_Bytes/Bytes_Impl_.h>
#include <_std/Type.h>
#include <apparentridges/Util.h>
#include <apparentridges/_ApparentRidges/Vec3_Impl_.h>
#include <apparentridges/_ApparentRidges/Face_Impl_.h>
#include <apparentridges/OBJParser.h>
#include <apparentridges/_ApparentRidges/Polyline_Impl_.h>
#include <apparentridges/PixelMap.h>
#include <apparentridges/SVGWriter.h>
#include <haxe/NativeStackTrace.h>
#include <haxe/ds/ArraySort.h>
#include <haxe/io/Input.h>
#include <hl/_NativeArray/NativeArray_Impl_.h>
#include <hl/_Type/Type_Impl_.h>
#include <hl/types/ArrayDyn.h>
#include <hl/types/_BytesMap/BytesMap_Impl_.h>
#include <sys/io/File.h>
#include <sys/io/FileInput.h>
extern hl_type t$String;
extern vbyte string$c6abe5e[];

// Globals
hl__$BaseType g$_hl_BaseType = 0;
hl__Class g$_hl_Class = 0;
$String g$_String = 0;
$Date g$_Date = 0;
$Std g$_Std = 0;
$Math g$_Math = 0;
String s$Can_t_add_ = 0;
String s$84c4047 = 0;
String s$_and_ = 0;
String s$9371d7a = 0;
String s$Invalid_unicode_char_ = 0;
String s$null = 0;
String s$ = 0;
hl__types__$ArrayAccess g$_hl_types_ArrayAccess = 0;
hl__types__$ArrayBase g$_hl_types_ArrayBase = 0;
hl__types__$ArrayBytes_hl_F32 g$_hl_types_ArrayBytes_hl_F32 = 0;
$StringBuf g$_StringBuf = 0;
$SysError g$_SysError = 0;
String s$SysError_ = 0;
$Sys g$_Sys = 0;
String s$68b329d = 0;
hl_bytes_map* g$__types__ = 0;
hl__$Enum g$_hl_Enum = 0;
apparentridges__$Ridge g$_apparentridges_Ridge = 0;
apparentridges__$BSphere g$_apparentridges_BSphere = 0;
apparentridges__$BBox g$_apparentridges_BBox = 0;
apparentridges__$BVHNode g$_apparentridges_BVHNode = 0;
apparentridges__$BVHTree g$_apparentridges_BVHTree = 0;
apparentridges__$Mesh g$_apparentridges_Mesh = 0;
haxe__$Log g$_haxe_Log = 0;
String s$computing_normals_ = 0;
String s$apparentridges_ApparentRidges_hx = 0;
String s$apparentridges_Mesh = 0;
String s$precompute = 0;
String s$computing_point_areas_ = 0;
String s$computing_adjacent_faces_ = 0;
String s$computing_curvatures_ = 0;
String s$computing_bounding_sphere_ = 0;
String s$computing_feature_size_ = 0;
String s$f2182c4 = 0;
String s$pre_computation_finished_ = 0;
apparentridges__$RayHit g$_apparentridges_RayHit = 0;
apparentridges__$Ray g$_apparentridges_Ray = 0;
apparentridges__$BVHPartition g$_apparentridges_BVHPartition = 0;
apparentridges__$BVHBucket g$_apparentridges_BVHBucket = 0;
String s$01abfc7 = 0;
String s$7215ee9 = 0;
String s$v = 0;
String s$f = 0;
String s$6666cd7 = 0;
apparentridges__$Line g$_apparentridges_Line = 0;
apparentridges__$Render g$_apparentridges_Render = 0;
String s$precomputing_mesh_properties_ = 0;
String s$apparentridges_Render = 0;
String s$apparentRidges = 0;
String s$generating_apparent_ridges_ = 0;
String s$cb2752f = 0;
String s$65d3534 = 0;
String s$de56809 = 0;
String s$buildPolylines = 0;
String s$polylines_built_ = 0;
String s$P3_ = 0;
String s$_255_ = 0;
String s$65e0ed6 = 0;
String s$_height_ = 0;
String s$50ed75a = 0;
String s$abc9c01 = 0;
String s$_y1_ = 0;
String s$_x2_ = 0;
String s$_y2_ = 0;
String s$_fill_none_stroke_rgb_ = 0;
String s$c0cb5f0 = 0;
String s$2411c92 = 0;
String s$_svg_ = 0;
String s$black = 0;
String s$_polyline_points_ = 0;
String s$rgb_ = 0;
String s$_fill_none_stroke_ = 0;
String s$c6abe5e = 0;
String s$_circle_cx_ = 0;
String s$_cy_ = 0;
String s$6c3ec5f = 0;
String s$6abd38a = 0;
String s$74b2e85 = 0;
String s$51ffe49 = 0;
String s$_stroke_rgb_ = 0;
String s$_g_ = 0;
haxe__$Exception g$_haxe_Exception = 0;
haxe__$ValueException g$_haxe_ValueException = 0;
String s$853ae90 = 0;
String s$fc763cb = 0;
String s$e265492 = 0;
String s$stack = 0;
String s$NativeStackTrace_callStack = 0;
haxe__io__$Encoding g$haxe_io_Encoding = 0;
haxe__io__$Bytes g$_haxe_io_Bytes = 0;
venum* g$haxe_io_Encoding_UTF8 = 0;
haxe__io__$BytesDataImpl g$_haxe_io_BytesDataImpl = 0;
haxe__io__$Eof g$_haxe_io_Eof = 0;
String s$Eof = 0;
haxe__io__$Output g$_haxe_io_Output = 0;
String s$Not_implemented = 0;
haxe__io__$Error g$haxe_io_Error = 0;
venum* g$haxe_io_Error_OutsideBounds = 0;
haxe__iterators__$ArrayIterator g$_haxe_iterators_ArrayIterator = 0;
haxe__iterators__$ArrayKeyValueIterator g$19142ef = 0;
hl__$NativeArrayIterator_Dynamic g$_hl_NativeArrayIterator_Dynamic = 0;
hl__$NativeArrayIterator_Int g$_hl_NativeArrayIterator_Int = 0;
hl__types__$BytesIterator_Float g$_hl_types_BytesIterator_Float = 0;
hl__types__$BytesKeyValueIterator_Float g$1e08d4d = 0;
String s$Invalid_array_index_ = 0;
hl__types__$BytesIterator_Int g$_hl_types_BytesIterator_Int = 0;
hl__types__$BytesKeyValueIterator_Int g$d250582 = 0;
hl__types__$BytesIterator_hl_F32 g$_hl_types_BytesIterator_hl_F32 = 0;
hl__types__$BytesKeyValueIterator_hl_F32 g$d4ac488 = 0;
hl__types__$BytesIterator_hl_UI16 g$_hl_types_BytesIterator_hl_UI16 = 0;
hl__types__$BytesKeyValueIterator_hl_UI16 g$4aef89b = 0;
hl__types__$ArrayDynIterator g$_hl_types_ArrayDynIterator = 0;
hl__types__$ArrayDynKeyValueIterator g$d29ca94 = 0;
hl__types__$ArrayObjIterator g$_hl_types_ArrayObjIterator = 0;
hl__types__$ArrayObjKeyValueIterator g$1813e5a = 0;
String s$2f43b42 = 0;
String s$Can_t_read_ = 0;
sys__io__$FileOutput g$_sys_io_FileOutput = 0;
String s$Can_t_open_ = 0;
String s$_for_writing = 0;
hl__CoreType g$_Float = 0;
String s$Float = 0;
hl__CoreType g$_Int = 0;
String s$Int = 0;
hl__CoreEnum g$_Bool = 0;
String s$Bool = 0;
hl__CoreType g$_Dynamic = 0;
String s$Dynamic = 0;
$StringTools g$_StringTools = 0;
hl___Bytes__$Bytes_Impl_ g$_hl__Bytes_Bytes_Impl_ = 0;
$Type g$_Type = 0;
apparentridges__$Util g$_apparentridges_Util = 0;
apparentridges___ApparentRidges__$Vec3_Impl_ g$8f3968b = 0;
apparentridges___ApparentRidges__$Face_Impl_ g$50a0721 = 0;
apparentridges__$OBJParser g$_apparentridges_OBJParser = 0;
apparentridges___ApparentRidges__$Polyline_Impl_ g$6121f78 = 0;
apparentridges__$PixelMap g$_apparentridges_PixelMap = 0;
apparentridges__$SVGWriter g$_apparentridges_SVGWriter = 0;
haxe__$NativeStackTrace g$_haxe_NativeStackTrace = 0;
haxe__ds__$ArraySort g$_haxe_ds_ArraySort = 0;
venum* g$haxe_io_Encoding_RawNative = 0;
venum* g$haxe_io_Error_Blocked = 0;
venum* g$haxe_io_Error_Overflow = 0;
haxe__io__$Input g$_haxe_io_Input = 0;
hl___NativeArray__$NativeArray_Impl_ g$7abf311 = 0;
hl___Type__$Type_Impl_ g$_hl__Type_Type_Impl_ = 0;
String s$Array = 0;
hl__types__$ArrayDyn g$_hl_types_ArrayDyn = 0;
String s$hl_types_ArrayDyn = 0;
hl__types___BytesMap__$BytesMap_Impl_ g$2c4fafe = 0;
sys__io__$File g$_sys_io_File = 0;
sys__io__$FileInput g$_sys_io_FileInput = 0;
static struct _String const_s$Can_t_add_ = {&t$String,(vbyte*)USTR("Can't add "),10};
static struct _String const_s$84c4047 = {&t$String,(vbyte*)USTR("("),1};
static struct _String const_s$_and_ = {&t$String,(vbyte*)USTR(") and "),6};
static struct _String const_s$9371d7a = {&t$String,(vbyte*)USTR(")"),1};
static struct _String const_s$Invalid_unicode_char_ = {&t$String,(vbyte*)USTR("Invalid unicode char "),21};
static struct _String const_s$null = {&t$String,(vbyte*)USTR("null"),4};
static struct _String const_s$ = {&t$String,(vbyte*)USTR(""),0};
static struct _String const_s$SysError_ = {&t$String,(vbyte*)USTR("SysError("),9};
static struct _String const_s$68b329d = {&t$String,(vbyte*)USTR("\n"),1};
static struct _String const_s$computing_normals_ = {&t$String,(vbyte*)USTR("computing normals..."),20};
static struct _String const_s$apparentridges_ApparentRidges_hx = {&t$String,(vbyte*)USTR("apparentridges/ApparentRidges.hx"),32};
static struct _String const_s$apparentridges_Mesh = {&t$String,(vbyte*)USTR("apparentridges.Mesh"),19};
static struct _String const_s$precompute = {&t$String,(vbyte*)USTR("precompute"),10};
static struct _String const_s$computing_point_areas_ = {&t$String,(vbyte*)USTR("computing point areas..."),24};
static struct _String const_s$computing_adjacent_faces_ = {&t$String,(vbyte*)USTR("computing adjacent faces..."),27};
static struct _String const_s$computing_curvatures_ = {&t$String,(vbyte*)USTR("computing curvatures..."),23};
static struct _String const_s$computing_bounding_sphere_ = {&t$String,(vbyte*)USTR("computing bounding sphere..."),28};
static struct _String const_s$computing_feature_size_ = {&t$String,(vbyte*)USTR("computing feature size..."),25};
static struct _String const_s$f2182c4 = {&t$String,(vbyte*)USTR("computing bounding volume hierarchy..."),38};
static struct _String const_s$pre_computation_finished_ = {&t$String,(vbyte*)USTR("pre-computation finished."),25};
static struct _String const_s$01abfc7 = {&t$String,(vbyte*)USTR("#"),1};
static struct _String const_s$7215ee9 = {&t$String,(vbyte*)USTR(" "),1};
static struct _String const_s$v = {&t$String,(vbyte*)USTR("v"),1};
static struct _String const_s$f = {&t$String,(vbyte*)USTR("f"),1};
static struct _String const_s$6666cd7 = {&t$String,(vbyte*)USTR("/"),1};
static struct _String const_s$precomputing_mesh_properties_ = {&t$String,(vbyte*)USTR("precomputing mesh properties..."),31};
static struct _String const_s$apparentridges_Render = {&t$String,(vbyte*)USTR("apparentridges.Render"),21};
static struct _String const_s$apparentRidges = {&t$String,(vbyte*)USTR("apparentRidges"),14};
static struct _String const_s$generating_apparent_ridges_ = {&t$String,(vbyte*)USTR("generating apparent ridges..."),29};
static struct _String const_s$cb2752f = {&t$String,(vbyte*)USTR("projecting apparent ridges onto 2D plane..."),43};
static struct _String const_s$65d3534 = {&t$String,(vbyte*)USTR("apparent ridges computation finished."),37};
static struct _String const_s$de56809 = {&t$String,(vbyte*)USTR("building polylines from ridge segments..."),41};
static struct _String const_s$buildPolylines = {&t$String,(vbyte*)USTR("buildPolylines"),14};
static struct _String const_s$polylines_built_ = {&t$String,(vbyte*)USTR("polylines built."),16};
static struct _String const_s$P3_ = {&t$String,(vbyte*)USTR("P3\n"),3};
static struct _String const_s$_255_ = {&t$String,(vbyte*)USTR("\n255\n"),5};
static struct _String const_s$65e0ed6 = {&t$String,(vbyte*)USTR("<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\""),61};
static struct _String const_s$_height_ = {&t$String,(vbyte*)USTR("\" height=\""),10};
static struct _String const_s$50ed75a = {&t$String,(vbyte*)USTR("\">\n"),3};
static struct _String const_s$abc9c01 = {&t$String,(vbyte*)USTR("  <line x1=\""),12};
static struct _String const_s$_y1_ = {&t$String,(vbyte*)USTR("\" y1=\""),6};
static struct _String const_s$_x2_ = {&t$String,(vbyte*)USTR("\" x2=\""),6};
static struct _String const_s$_y2_ = {&t$String,(vbyte*)USTR("\" y2=\""),6};
static struct _String const_s$_fill_none_stroke_rgb_ = {&t$String,(vbyte*)USTR("\" fill=\"none\" stroke=\"rgb("),26};
static struct _String const_s$c0cb5f0 = {&t$String,(vbyte*)USTR(","),1};
static struct _String const_s$2411c92 = {&t$String,(vbyte*)USTR(")\" stroke-width=\"1\" stroke-linecap=\"round\"/>\n"),45};
static struct _String const_s$_svg_ = {&t$String,(vbyte*)USTR("</svg>\n"),7};
static struct _String const_s$black = {&t$String,(vbyte*)USTR("black"),5};
static struct _String const_s$_polyline_points_ = {&t$String,(vbyte*)USTR("  <polyline points=\""),20};
static struct _String const_s$rgb_ = {&t$String,(vbyte*)USTR("rgb("),4};
static struct _String const_s$_fill_none_stroke_ = {&t$String,(vbyte*)USTR("\" fill=\"none\" stroke=\""),22};
static struct _String const_s$c6abe5e = {&t$String,(vbyte*)string$c6abe5e,68};
static struct _String const_s$_circle_cx_ = {&t$String,(vbyte*)USTR("  <circle cx=\""),14};
static struct _String const_s$_cy_ = {&t$String,(vbyte*)USTR("\" cy=\""),6};
static struct _String const_s$6c3ec5f = {&t$String,(vbyte*)USTR("\" r=\"0.5\" stroke=\"none\" fill=\"rgb("),34};
static struct _String const_s$6abd38a = {&t$String,(vbyte*)USTR(")\"/>\n"),5};
static struct _String const_s$74b2e85 = {&t$String,(vbyte*)USTR("  <g fill=\"none\" stroke-width=\"1\">\n"),35};
static struct _String const_s$51ffe49 = {&t$String,(vbyte*)USTR("    <line x1=\""),14};
static struct _String const_s$_stroke_rgb_ = {&t$String,(vbyte*)USTR("\" stroke=\"rgb("),14};
static struct _String const_s$_g_ = {&t$String,(vbyte*)USTR("  </g>\n"),7};
static struct _String const_s$853ae90 = {&t$String,(vbyte*)USTR(":"),1};
static struct _String const_s$fc763cb = {&t$String,(vbyte*)USTR(", "),2};
static struct _String const_s$e265492 = {&t$String,(vbyte*)USTR(": "),2};
static struct _String const_s$stack = {&t$String,(vbyte*)USTR("stack"),5};
static struct _String const_s$NativeStackTrace_callStack = {&t$String,(vbyte*)USTR("NativeStackTrace.callStack"),26};
static struct _String const_s$Eof = {&t$String,(vbyte*)USTR("Eof"),3};
static struct _String const_s$Not_implemented = {&t$String,(vbyte*)USTR("Not implemented"),15};
static struct _String const_s$Invalid_array_index_ = {&t$String,(vbyte*)USTR("Invalid array index "),20};
static struct _String const_s$2f43b42 = {&t$String,(vbyte*)USTR("..."),3};
static struct _String const_s$Can_t_read_ = {&t$String,(vbyte*)USTR("Can't read "),11};
static struct _String const_s$Can_t_open_ = {&t$String,(vbyte*)USTR("Can't open "),11};
static struct _String const_s$_for_writing = {&t$String,(vbyte*)USTR(" for writing"),12};
static struct _String const_s$Float = {&t$String,(vbyte*)USTR("Float"),5};
static struct _String const_s$Int = {&t$String,(vbyte*)USTR("Int"),3};
static struct _String const_s$Bool = {&t$String,(vbyte*)USTR("Bool"),4};
static struct _String const_s$Dynamic = {&t$String,(vbyte*)USTR("Dynamic"),7};
static struct _String const_s$Array = {&t$String,(vbyte*)USTR("Array"),5};
static struct _String const_s$hl_types_ArrayDyn = {&t$String,(vbyte*)USTR("hl.types.ArrayDyn"),17};

void hl_init_roots() {
	s$Can_t_add_ = &const_s$Can_t_add_;
	s$84c4047 = &const_s$84c4047;
	s$_and_ = &const_s$_and_;
	s$9371d7a = &const_s$9371d7a;
	s$Invalid_unicode_char_ = &const_s$Invalid_unicode_char_;
	s$null = &const_s$null;
	s$ = &const_s$;
	s$SysError_ = &const_s$SysError_;
	s$68b329d = &const_s$68b329d;
	s$computing_normals_ = &const_s$computing_normals_;
	s$apparentridges_ApparentRidges_hx = &const_s$apparentridges_ApparentRidges_hx;
	s$apparentridges_Mesh = &const_s$apparentridges_Mesh;
	s$precompute = &const_s$precompute;
	s$computing_point_areas_ = &const_s$computing_point_areas_;
	s$computing_adjacent_faces_ = &const_s$computing_adjacent_faces_;
	s$computing_curvatures_ = &const_s$computing_curvatures_;
	s$computing_bounding_sphere_ = &const_s$computing_bounding_sphere_;
	s$computing_feature_size_ = &const_s$computing_feature_size_;
	s$f2182c4 = &const_s$f2182c4;
	s$pre_computation_finished_ = &const_s$pre_computation_finished_;
	s$01abfc7 = &const_s$01abfc7;
	s$7215ee9 = &const_s$7215ee9;
	s$v = &const_s$v;
	s$f = &const_s$f;
	s$6666cd7 = &const_s$6666cd7;
	s$precomputing_mesh_properties_ = &const_s$precomputing_mesh_properties_;
	s$apparentridges_Render = &const_s$apparentridges_Render;
	s$apparentRidges = &const_s$apparentRidges;
	s$generating_apparent_ridges_ = &const_s$generating_apparent_ridges_;
	s$cb2752f = &const_s$cb2752f;
	s$65d3534 = &const_s$65d3534;
	s$de56809 = &const_s$de56809;
	s$buildPolylines = &const_s$buildPolylines;
	s$polylines_built_ = &const_s$polylines_built_;
	s$P3_ = &const_s$P3_;
	s$_255_ = &const_s$_255_;
	s$65e0ed6 = &const_s$65e0ed6;
	s$_height_ = &const_s$_height_;
	s$50ed75a = &const_s$50ed75a;
	s$abc9c01 = &const_s$abc9c01;
	s$_y1_ = &const_s$_y1_;
	s$_x2_ = &const_s$_x2_;
	s$_y2_ = &const_s$_y2_;
	s$_fill_none_stroke_rgb_ = &const_s$_fill_none_stroke_rgb_;
	s$c0cb5f0 = &const_s$c0cb5f0;
	s$2411c92 = &const_s$2411c92;
	s$_svg_ = &const_s$_svg_;
	s$black = &const_s$black;
	s$_polyline_points_ = &const_s$_polyline_points_;
	s$rgb_ = &const_s$rgb_;
	s$_fill_none_stroke_ = &const_s$_fill_none_stroke_;
	s$c6abe5e = &const_s$c6abe5e;
	s$_circle_cx_ = &const_s$_circle_cx_;
	s$_cy_ = &const_s$_cy_;
	s$6c3ec5f = &const_s$6c3ec5f;
	s$6abd38a = &const_s$6abd38a;
	s$74b2e85 = &const_s$74b2e85;
	s$51ffe49 = &const_s$51ffe49;
	s$_stroke_rgb_ = &const_s$_stroke_rgb_;
	s$_g_ = &const_s$_g_;
	s$853ae90 = &const_s$853ae90;
	s$fc763cb = &const_s$fc763cb;
	s$e265492 = &const_s$e265492;
	s$stack = &const_s$stack;
	s$NativeStackTrace_callStack = &const_s$NativeStackTrace_callStack;
	s$Eof = &const_s$Eof;
	s$Not_implemented = &const_s$Not_implemented;
	s$Invalid_array_index_ = &const_s$Invalid_array_index_;
	s$2f43b42 = &const_s$2f43b42;
	s$Can_t_read_ = &const_s$Can_t_read_;
	s$Can_t_open_ = &const_s$Can_t_open_;
	s$_for_writing = &const_s$_for_writing;
	s$Float = &const_s$Float;
	s$Int = &const_s$Int;
	s$Bool = &const_s$Bool;
	s$Dynamic = &const_s$Dynamic;
	s$Array = &const_s$Array;
	s$hl_types_ArrayDyn = &const_s$hl_types_ArrayDyn;
	hl_add_root((void**)&g$_hl_BaseType);
	hl_add_root((void**)&g$_hl_Class);
	hl_add_root((void**)&g$_String);
	hl_add_root((void**)&g$_Date);
	hl_add_root((void**)&g$_Std);
	hl_add_root((void**)&g$_Math);
	hl_add_root((void**)&g$_hl_types_ArrayAccess);
	hl_add_root((void**)&g$_hl_types_ArrayBase);
	hl_add_root((void**)&g$_hl_types_ArrayBytes_hl_F32);
	hl_add_root((void**)&g$_StringBuf);
	hl_add_root((void**)&g$_SysError);
	hl_add_root((void**)&g$_Sys);
	hl_add_root((void**)&g$__types__);
	hl_add_root((void**)&g$_hl_Enum);
	hl_add_root((void**)&g$_apparentridges_Ridge);
	hl_add_root((void**)&g$_apparentridges_BSphere);
	hl_add_root((void**)&g$_apparentridges_BBox);
	hl_add_root((void**)&g$_apparentridges_BVHNode);
	hl_add_root((void**)&g$_apparentridges_BVHTree);
	hl_add_root((void**)&g$_apparentridges_Mesh);
	hl_add_root((void**)&g$_haxe_Log);
	hl_add_root((void**)&g$_apparentridges_RayHit);
	hl_add_root((void**)&g$_apparentridges_Ray);
	hl_add_root((void**)&g$_apparentridges_BVHPartition);
	hl_add_root((void**)&g$_apparentridges_BVHBucket);
	hl_add_root((void**)&g$_apparentridges_Line);
	hl_add_root((void**)&g$_apparentridges_Render);
	hl_add_root((void**)&g$_haxe_Exception);
	hl_add_root((void**)&g$_haxe_ValueException);
	hl_add_root((void**)&g$haxe_io_Encoding);
	hl_add_root((void**)&g$_haxe_io_Bytes);
	hl_add_root((void**)&g$haxe_io_Encoding_UTF8);
	hl_add_root((void**)&g$_haxe_io_BytesDataImpl);
	hl_add_root((void**)&g$_haxe_io_Eof);
	hl_add_root((void**)&g$_haxe_io_Output);
	hl_add_root((void**)&g$haxe_io_Error);
	hl_add_root((void**)&g$haxe_io_Error_OutsideBounds);
	hl_add_root((void**)&g$_haxe_iterators_ArrayIterator);
	hl_add_root((void**)&g$19142ef);
	hl_add_root((void**)&g$_hl_NativeArrayIterator_Dynamic);
	hl_add_root((void**)&g$_hl_NativeArrayIterator_Int);
	hl_add_root((void**)&g$_hl_types_BytesIterator_Float);
	hl_add_root((void**)&g$1e08d4d);
	hl_add_root((void**)&g$_hl_types_BytesIterator_Int);
	hl_add_root((void**)&g$d250582);
	hl_add_root((void**)&g$_hl_types_BytesIterator_hl_F32);
	hl_add_root((void**)&g$d4ac488);
	hl_add_root((void**)&g$_hl_types_BytesIterator_hl_UI16);
	hl_add_root((void**)&g$4aef89b);
	hl_add_root((void**)&g$_hl_types_ArrayDynIterator);
	hl_add_root((void**)&g$d29ca94);
	hl_add_root((void**)&g$_hl_types_ArrayObjIterator);
	hl_add_root((void**)&g$1813e5a);
	hl_add_root((void**)&g$_sys_io_FileOutput);
	hl_add_root((void**)&g$_Float);
	hl_add_root((void**)&g$_Int);
	hl_add_root((void**)&g$_Bool);
	hl_add_root((void**)&g$_Dynamic);
	hl_add_root((void**)&g$_StringTools);
	hl_add_root((void**)&g$_hl__Bytes_Bytes_Impl_);
	hl_add_root((void**)&g$_Type);
	hl_add_root((void**)&g$_apparentridges_Util);
	hl_add_root((void**)&g$8f3968b);
	hl_add_root((void**)&g$50a0721);
	hl_add_root((void**)&g$_apparentridges_OBJParser);
	hl_add_root((void**)&g$6121f78);
	hl_add_root((void**)&g$_apparentridges_PixelMap);
	hl_add_root((void**)&g$_apparentridges_SVGWriter);
	hl_add_root((void**)&g$_haxe_NativeStackTrace);
	hl_add_root((void**)&g$_haxe_ds_ArraySort);
	hl_add_root((void**)&g$haxe_io_Encoding_RawNative);
	hl_add_root((void**)&g$haxe_io_Error_Blocked);
	hl_add_root((void**)&g$haxe_io_Error_Overflow);
	hl_add_root((void**)&g$_haxe_io_Input);
	hl_add_root((void**)&g$7abf311);
	hl_add_root((void**)&g$_hl__Type_Type_Impl_);
	hl_add_root((void**)&g$_hl_types_ArrayDyn);
	hl_add_root((void**)&g$2c4fafe);
	hl_add_root((void**)&g$_sys_io_File);
	hl_add_root((void**)&g$_sys_io_FileInput);
}
// \" stroke-width=\"1\" stroke-linecap=\"round\" stroke-linejoin=\"r...
vbyte string$c6abe5e[] = {34,0,32,0,115,0,116,0,114,0,111,0,107,0,101,0,45,0,119,0,105,0,100,0,116,0,104,0,61,0,34,0,49,0,34,0,32,0,115,0,116,0,114,0,111,0,107,0,101,0,45,0,108,0,105,0,110,0,101,0,99,0,97,0,112,0,61,0,34,0,114,0,111,0,117,0,110,0,100,0,34,0,32,0,115,0,116,0,114,0,111,0,107,0,101,0,45,0,108,0,105,0,110,0,101,0,106,0,111,0,105,0,110,0,61,0,34,0,114,0,111,0,117,0,110,0,100\
	,0,34,0,47,0,62,0,10,0,0,0};
