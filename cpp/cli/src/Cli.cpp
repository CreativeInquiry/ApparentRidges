// Generated by Haxe 4.1.3
#include <hxcpp.h>

#ifndef INCLUDED_95f339a1d026d52c
#define INCLUDED_95f339a1d026d52c
#include "hxMath.h"
#endif
#ifndef INCLUDED_Cli
#include <Cli.h>
#endif
#ifndef INCLUDED_Std
#include <Std.h>
#endif
#ifndef INCLUDED_StringTools
#include <StringTools.h>
#endif
#ifndef INCLUDED_Sys
#include <Sys.h>
#endif
#ifndef INCLUDED_apparentridges_Mesh
#include <apparentridges/Mesh.h>
#endif
#ifndef INCLUDED_apparentridges_OBJParser
#include <apparentridges/OBJParser.h>
#endif
#ifndef INCLUDED_apparentridges_PixelMap
#include <apparentridges/PixelMap.h>
#endif
#ifndef INCLUDED_apparentridges_Render
#include <apparentridges/Render.h>
#endif
#ifndef INCLUDED_apparentridges_SVGWriter
#include <apparentridges/SVGWriter.h>
#endif
#ifndef INCLUDED_haxe_IMap
#include <haxe/IMap.h>
#endif
#ifndef INCLUDED_haxe_Log
#include <haxe/Log.h>
#endif
#ifndef INCLUDED_haxe_ds_StringMap
#include <haxe/ds/StringMap.h>
#endif
#ifndef INCLUDED_sys_io_File
#include <sys/io/File.h>
#endif

HX_LOCAL_STACK_FRAME(_hx_pos_0ede3bd59b9ae14d_5_main,"Cli","main",0x437fc1e7,"Cli.main","Cli.hx",5,0x3eba7d3e)
static const Float _hx_array_data_00333580_8[] = {
	(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,
};
static const Float _hx_array_data_00333580_9[] = {
	(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,(Float)0,(Float)0,(Float)0,(Float)0,(Float)1,
};

void Cli_obj::__construct() { }

Dynamic Cli_obj::__CreateEmpty() { return new Cli_obj; }

void *Cli_obj::_hx_vtable = 0;

Dynamic Cli_obj::__Create(::hx::DynamicArray inArgs)
{
	::hx::ObjectPtr< Cli_obj > _hx_result = new Cli_obj();
	_hx_result->__construct();
	return _hx_result;
}

bool Cli_obj::_hx_isInstanceOf(int inClassId) {
	return inClassId==(int)0x00000001 || inClassId==(int)0x7ed77a14;
}

int Cli_obj::main(){
            	HX_GC_STACKFRAME(&_hx_pos_0ede3bd59b9ae14d_5_main)
HXLINE(  16)		 ::haxe::ds::StringMap _g =  ::haxe::ds::StringMap_obj::__alloc( HX_CTX );
HXDLIN(  16)		_g->set(HX_("-o",a2,27,00,00),HX_("--output",61,fe,a4,69));
HXDLIN(  16)		_g->set(HX_("-w",aa,27,00,00),HX_("--width",a6,d2,b7,17));
HXDLIN(  16)		_g->set(HX_("-h",9b,27,00,00),HX_("--height",47,f7,6f,5f));
HXDLIN(  16)		 ::haxe::ds::StringMap abbrev = _g;
HXLINE(  17)		::Array< ::String > argv = ::Sys_obj::args();
HXLINE(  18)		 ::haxe::ds::StringMap _g1 =  ::haxe::ds::StringMap_obj::__alloc( HX_CTX );
HXDLIN(  18)		_g1->set(HX_("--width",a6,d2,b7,17),HX_("800",38,a8,2a,00));
HXDLIN(  18)		_g1->set(HX_("--height",47,f7,6f,5f),HX_("800",38,a8,2a,00));
HXDLIN(  18)		{
HXLINE(  21)			::String id = HX_("",00,00,00,00);
HXDLIN(  21)			{
HXLINE(  21)				id = (id + ::String::fromCharCode((::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(26) ))) + 97)));
HXDLIN(  21)				id = (id + ::String::fromCharCode((::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(26) ))) + 97)));
HXDLIN(  21)				id = (id + ::String::fromCharCode((::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(26) ))) + 97)));
HXDLIN(  21)				id = (id + ::String::fromCharCode((::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(26) ))) + 97)));
            			}
HXLINE(  18)			_g1->set(HX_("--output",61,fe,a4,69),((HX_("out-",df,db,b7,49) + id) + HX_(".svg",f6,7a,bf,1e)));
            		}
HXDLIN(  18)		_g1->set(HX_("--transform",0c,1a,c3,b9),HX_("auto()",30,fc,80,73));
HXDLIN(  18)		_g1->set(HX_("--cull",12,9d,e0,28),HX_("true",4e,a7,03,4d));
HXDLIN(  18)		_g1->set(HX_("--thresh",dc,79,dd,eb),HX_("0.1",73,94,24,00));
HXDLIN(  18)		_g1->set(HX_("--mode",c3,2a,78,2f),HX_("gradients",83,78,13,cd));
HXDLIN(  18)		_g1->set(HX_("--verbose",22,5c,07,94),HX_("1",31,00,00,00));
HXDLIN(  18)		 ::haxe::ds::StringMap kwargs = _g1;
HXLINE(  28)		 ::haxe::ds::StringMap _g2 =  ::haxe::ds::StringMap_obj::__alloc( HX_CTX );
HXDLIN(  28)		_g2->set(HX_("--width",a6,d2,b7,17),HX_("width of output image",8b,63,44,b7));
HXDLIN(  28)		_g2->set(HX_("--height",47,f7,6f,5f),HX_("height of output image",6c,47,eb,4e));
HXDLIN(  28)		_g2->set(HX_("--output",61,fe,a4,69),HX_("output filename\n(randomly generated if not specified)",4c,ba,7e,e2));
HXDLIN(  28)		_g2->set(HX_("--transform",0c,1a,c3,b9),HX_("place the model before rendering\n(multiple commands are left-multiplied in order)\navailable commands:\nfocal(100) translate(1,2,3) scale(4,4,4)\nrotateX(1) rotateY(2) rotateZ(3);\nmatrix(11,12,13,14,21,22,...,43,44)\nrotation angles are in degrees\nuse auto() or auto(zFactor,fFactor) for automatic\nplacement. (cam is fixed at (0,0,0) pointing at +Z,\nfocal() needs to be specified for manual placement)",e9,e5,f8,38));
HXDLIN(  28)		_g2->set(HX_("--cull",12,9d,e0,28),HX_("don't draw occluded faces\noptions: true/false/custom float value (e.g. 1.5)\nthe float being a multiplier to the 'epsilon',\nallowing ridges just barely occluded to show up",cb,42,dd,db));
HXDLIN(  28)		_g2->set(HX_("--thresh",dc,79,dd,eb),HX_("apparent ridge threshold\nsmaller => more detailed, larger => cleaner",49,04,74,a0));
HXDLIN(  28)		_g2->set(HX_("--mode",c3,2a,78,2f),HX_("visualization technique, options:\nvertices, edges, \nlines (disconnected ridge segments),\npolylines (connected ridges via post-proc),\ngradients (nice render showing ridge strength),\npixels (raster out: depth, norm., curv. maps)",34,dc,e1,df));
HXDLIN(  28)		_g2->set(HX_("--verbose",22,5c,07,94),HX_("verbosity: 0 => errors only, 1 => logs",e8,1c,f1,b3));
HXDLIN(  28)		 ::haxe::ds::StringMap desc = _g2;
HXLINE(  39)		bool _hx_tmp;
HXDLIN(  39)		if ((argv->length != 0)) {
HXLINE(  39)			_hx_tmp = (argv->__get(0) == HX_("--help",21,8a,22,2c));
            		}
            		else {
HXLINE(  39)			_hx_tmp = true;
            		}
HXDLIN(  39)		if (_hx_tmp) {
HXLINE(  40)			::Sys_obj::println(HX_("\n       ___           ___      \n      /\\  \\         /\\  \\     \n     /  \\  \\       /  \\  \\    \n    / /\\ \\  \\     / /\\ \\  \\   \n   /  \\~\\ \\  \\   /  \\~\\ \\  \\  \n  / /\\ \\ \\ \\__\\ / /\\ \\ \\ \\__\\ \n  \\/__\\ \\/ /  / \\/_|  \\/ /  / \n       \\  /  /     | |  /  /  \n       / /  /      | |\\/__/   \n      / /  /       | |  |     \n      \\/__/ PPARENT \\|__| IDGES (CLI)\n\nRender line drawings of 3D meshes.\n\nusage:   apparentridges [options] [.obj files...]\n\noption        default       description\n",aa,b5,b3,a8));
HXLINE(  59)			{
HXLINE(  59)				 ::Dynamic k = kwargs->keys();
HXDLIN(  59)				while(( (bool)(k->__Field(HX_("hasNext",6d,a5,46,18),::hx::paccDynamic)()) )){
HXLINE(  59)					::String k1 = ( (::String)(k->__Field(HX_("next",f3,84,02,49),::hx::paccDynamic)()) );
HXLINE(  60)					::String defau = kwargs->get_string(k1);
HXLINE(  61)					::String ln = k1;
HXLINE(  62)					{
HXLINE(  62)						 ::Dynamic a = abbrev->keys();
HXDLIN(  62)						while(( (bool)(a->__Field(HX_("hasNext",6d,a5,46,18),::hx::paccDynamic)()) )){
HXLINE(  62)							::String a1 = ( (::String)(a->__Field(HX_("next",f3,84,02,49),::hx::paccDynamic)()) );
HXLINE(  63)							if ((abbrev->get_string(a1) == k1)) {
HXLINE(  64)								ln = (ln + (HX_(", ",74,26,00,00) + a1));
            							}
            						}
            					}
HXLINE(  67)					ln = ::StringTools_obj::rpad(ln,HX_(" ",20,00,00,00),14);
HXLINE(  68)					ln = (ln + defau);
HXLINE(  69)					ln = ::StringTools_obj::rpad(ln,HX_(" ",20,00,00,00),28);
HXLINE(  70)					if (desc->exists(k1)) {
HXLINE(  71)						::String ln1 = desc->get_string(k1);
HXDLIN(  71)						::Array< ::String > _g = ::Array_obj< ::String >::__new(0);
HXDLIN(  71)						{
HXLINE(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
HXDLIN(  71)							_g->push(HX_(" ",20,00,00,00));
            						}
HXDLIN(  71)						ln = (ln + ::StringTools_obj::replace(ln1,HX_("\n",0a,00,00,00),(HX_("\n",0a,00,00,00) + _g->join(HX_("",00,00,00,00)))));
            					}
HXLINE(  73)					ln = (ln + HX_("\n",0a,00,00,00));
HXLINE(  74)					::Sys_obj::println(ln);
            				}
            			}
HXLINE(  77)			return 0;
            		}
HXLINE(  80)		::Array< ::String > args = ::Array_obj< ::String >::__new(0);
HXLINE(  82)		int i = 0;
HXLINE(  83)		while((i < argv->length)){
HXLINE(  84)			if ((argv->__get(i).length == 0)) {
HXLINE(  85)				i = (i + 1);
HXLINE(  86)				continue;
            			}
HXLINE(  88)			if ((argv->__get(i).charAt(0) == HX_("-",2d,00,00,00))) {
HXLINE(  89)				if (((i + 1) >= argv->length)) {
HXLINE(  90)					::haxe::Log_obj::trace((HX_("[error] option without value: ",c4,c3,d8,d2) + argv->__get(i)),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),90,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
HXLINE(  91)					return 1;
            				}
HXLINE(  93)				::String key = argv->__get(i);
HXLINE(  94)				::String val = argv->__get((i + 1));
HXLINE(  95)				if (abbrev->exists(key)) {
HXLINE(  96)					key = abbrev->get_string(key);
            				}
HXLINE(  98)				{
HXLINE(  98)					::String v = ::StringTools_obj::replace(val,HX_(" ",20,00,00,00),HX_("",00,00,00,00));
HXDLIN(  98)					kwargs->set(key,v);
            				}
HXLINE(  99)				i = (i + 2);
HXLINE( 100)				continue;
            			}
HXLINE( 102)			args->push(argv->__get(i));
HXLINE( 103)			i = (i + 1);
            		}
HXLINE( 105)		if ((args->length < 1)) {
HXLINE( 106)			::haxe::Log_obj::trace(HX_("[warn] no input file.",c1,e9,63,45),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),106,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
HXLINE( 107)			return 0;
            		}
HXLINE( 117)		int verb = ( (int)(::Std_obj::parseInt(kwargs->get_string(HX_("--verbose",22,5c,07,94)))) );
HXLINE( 119)		{
HXLINE( 119)			int _g3 = 0;
HXDLIN( 119)			while((_g3 < args->length)){
HXLINE( 119)				::String inp = args->__get(_g3);
HXDLIN( 119)				_g3 = (_g3 + 1);
HXLINE( 120)				if ((verb > 0)) {
HXLINE( 120)					::haxe::Log_obj::trace(((HX_("[info] working on ",06,ae,11,b7) + inp) + HX_("...",ee,0f,23,00)),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),120,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            				}
HXLINE( 121)				 ::apparentridges::Mesh mesh = ::apparentridges::OBJParser_obj::fromFile(inp);
HXLINE( 123)				 ::Dynamic render = ::Std_obj::parseInt(kwargs->get_string(HX_("--width",a6,d2,b7,17)));
HXLINE( 122)				 ::apparentridges::Render render1 =  ::apparentridges::Render_obj::__alloc( HX_CTX ,mesh,( (int)(render) ),( (int)(::Std_obj::parseInt(kwargs->get_string(HX_("--height",47,f7,6f,5f)))) ));
HXLINE( 125)				render1->setVerbose(verb);
HXLINE( 127)				::Array< Float > m = ::Array_obj< Float >::fromData( _hx_array_data_00333580_8,16);
HXLINE( 128)				{
HXLINE( 128)					int _g = 0;
HXDLIN( 128)					::Array< ::String > _g1 = kwargs->get_string(HX_("--transform",0c,1a,c3,b9)).split(HX_(")",29,00,00,00));
HXDLIN( 128)					while((_g < _g1->length)){
HXLINE( 128)						::String t = _g1->__get(_g);
HXDLIN( 128)						_g = (_g + 1);
HXLINE( 129)						if ((t.length == 0)) {
HXLINE( 130)							continue;
            						}
HXLINE( 132)						::Array< Float > _g2 = ::Array_obj< Float >::__new(0);
HXDLIN( 132)						{
HXLINE( 132)							int _g3 = 0;
HXDLIN( 132)							::Array< ::String > _g4 = t.split(HX_("(",28,00,00,00))->__get(1).split(HX_(")",29,00,00,00))->__get(0).split(HX_(",",2c,00,00,00));
HXDLIN( 132)							while((_g3 < _g4->length)){
HXLINE( 132)								::String x = _g4->__get(_g3);
HXDLIN( 132)								_g3 = (_g3 + 1);
HXDLIN( 132)								_g2->push(::Std_obj::parseFloat(x));
            							}
            						}
HXDLIN( 132)						::Array< Float > _g5 = ::Array_obj< Float >::__new(0);
HXDLIN( 132)						{
HXLINE( 132)							int _g6 = 0;
HXDLIN( 132)							::Array< Float > _g7 = _g2;
HXDLIN( 132)							while((_g6 < _g7->length)){
HXLINE( 132)								Float v = _g7->__get(_g6);
HXDLIN( 132)								_g6 = (_g6 + 1);
HXDLIN( 132)								if (!(::Math_obj::isNaN(v))) {
HXLINE( 132)									_g5->push(v);
            								}
            							}
            						}
HXDLIN( 132)						::Array< Float > nums = _g5;
HXLINE( 134)						if (::StringTools_obj::startsWith(t,HX_("rotateX",9d,49,1d,f1))) {
HXLINE( 135)							Float a = ((nums->__get(0) * ::Math_obj::PI) / ( (Float)(180) ));
HXDLIN( 135)							Float A_0 = ( (Float)(1) );
HXDLIN( 135)							Float A_1 = ( (Float)(0) );
HXDLIN( 135)							Float A_2 = ( (Float)(0) );
HXDLIN( 135)							Float A_3 = ( (Float)(0) );
HXDLIN( 135)							Float A_4 = ( (Float)(0) );
HXDLIN( 135)							Float A_5 = ::Math_obj::cos(a);
HXDLIN( 135)							Float A_6 = -(::Math_obj::sin(a));
HXDLIN( 135)							Float A_7 = ( (Float)(0) );
HXDLIN( 135)							Float A_8 = ( (Float)(0) );
HXDLIN( 135)							Float A_9 = ::Math_obj::sin(a);
HXDLIN( 135)							Float A_10 = ::Math_obj::cos(a);
HXDLIN( 135)							Float A_11 = ( (Float)(0) );
HXDLIN( 135)							Float A_12 = ( (Float)(0) );
HXDLIN( 135)							Float A_13 = ( (Float)(0) );
HXDLIN( 135)							Float A_14 = ( (Float)(0) );
HXDLIN( 135)							Float A_15 = ( (Float)(1) );
HXDLIN( 135)							m = ::Array_obj< Float >::__new(16)->init(0,((((A_0 * m->__get(0)) + (A_1 * m->__get(4))) + (A_2 * m->__get(8))) + (A_3 * m->__get(12))))->init(1,((((A_0 * m->__get(1)) + (A_1 * m->__get(5))) + (A_2 * m->__get(9))) + (A_3 * m->__get(13))))->init(2,((((A_0 * m->__get(2)) + (A_1 * m->__get(6))) + (A_2 * m->__get(10))) + (A_3 * m->__get(14))))->init(3,((((A_0 * m->__get(3)) + (A_1 * m->__get(7))) + (A_2 * m->__get(11))) + (A_3 * m->__get(15))))->init(4,((((A_4 * m->__get(0)) + (A_5 * m->__get(4))) + (A_6 * m->__get(8))) + (A_7 * m->__get(12))))->init(5,((((A_4 * m->__get(1)) + (A_5 * m->__get(5))) + (A_6 * m->__get(9))) + (A_7 * m->__get(13))))->init(6,((((A_4 * m->__get(2)) + (A_5 * m->__get(6))) + (A_6 * m->__get(10))) + (A_7 * m->__get(14))))->init(7,((((A_4 * m->__get(3)) + (A_5 * m->__get(7))) + (A_6 * m->__get(11))) + (A_7 * m->__get(15))))->init(8,((((A_8 * m->__get(0)) + (A_9 * m->__get(4))) + (A_10 * m->__get(8))) + (A_11 * m->__get(12))))->init(9,((((A_8 * m->__get(1)) + (A_9 * m->__get(5))) + (A_10 * m->__get(9))) + (A_11 * m->__get(13))))->init(10,((((A_8 * m->__get(2)) + (A_9 * m->__get(6))) + (A_10 * m->__get(10))) + (A_11 * m->__get(14))))->init(11,((((A_8 * m->__get(3)) + (A_9 * m->__get(7))) + (A_10 * m->__get(11))) + (A_11 * m->__get(15))))->init(12,((((A_12 * m->__get(0)) + (A_13 * m->__get(4))) + (A_14 * m->__get(8))) + (A_15 * m->__get(12))))->init(13,((((A_12 * m->__get(1)) + (A_13 * m->__get(5))) + (A_14 * m->__get(9))) + (A_15 * m->__get(13))))->init(14,((((A_12 * m->__get(2)) + (A_13 * m->__get(6))) + (A_14 * m->__get(10))) + (A_15 * m->__get(14))))->init(15,((((A_12 * m->__get(3)) + (A_13 * m->__get(7))) + (A_14 * m->__get(11))) + (A_15 * m->__get(15))));
            						}
            						else {
HXLINE( 136)							if (::StringTools_obj::startsWith(t,HX_("rotateY",9e,49,1d,f1))) {
HXLINE( 137)								Float a = ((nums->__get(0) * ::Math_obj::PI) / ( (Float)(180) ));
HXDLIN( 137)								Float A_0 = ::Math_obj::cos(a);
HXDLIN( 137)								Float A_1 = ( (Float)(0) );
HXDLIN( 137)								Float A_2 = ::Math_obj::sin(a);
HXDLIN( 137)								Float A_3 = ( (Float)(0) );
HXDLIN( 137)								Float A_4 = ( (Float)(0) );
HXDLIN( 137)								Float A_5 = ( (Float)(1) );
HXDLIN( 137)								Float A_6 = ( (Float)(0) );
HXDLIN( 137)								Float A_7 = ( (Float)(0) );
HXDLIN( 137)								Float A_8 = -(::Math_obj::sin(a));
HXDLIN( 137)								Float A_9 = ( (Float)(0) );
HXDLIN( 137)								Float A_10 = ::Math_obj::cos(a);
HXDLIN( 137)								Float A_11 = ( (Float)(0) );
HXDLIN( 137)								Float A_12 = ( (Float)(0) );
HXDLIN( 137)								Float A_13 = ( (Float)(0) );
HXDLIN( 137)								Float A_14 = ( (Float)(0) );
HXDLIN( 137)								Float A_15 = ( (Float)(1) );
HXDLIN( 137)								m = ::Array_obj< Float >::__new(16)->init(0,((((A_0 * m->__get(0)) + (A_1 * m->__get(4))) + (A_2 * m->__get(8))) + (A_3 * m->__get(12))))->init(1,((((A_0 * m->__get(1)) + (A_1 * m->__get(5))) + (A_2 * m->__get(9))) + (A_3 * m->__get(13))))->init(2,((((A_0 * m->__get(2)) + (A_1 * m->__get(6))) + (A_2 * m->__get(10))) + (A_3 * m->__get(14))))->init(3,((((A_0 * m->__get(3)) + (A_1 * m->__get(7))) + (A_2 * m->__get(11))) + (A_3 * m->__get(15))))->init(4,((((A_4 * m->__get(0)) + (A_5 * m->__get(4))) + (A_6 * m->__get(8))) + (A_7 * m->__get(12))))->init(5,((((A_4 * m->__get(1)) + (A_5 * m->__get(5))) + (A_6 * m->__get(9))) + (A_7 * m->__get(13))))->init(6,((((A_4 * m->__get(2)) + (A_5 * m->__get(6))) + (A_6 * m->__get(10))) + (A_7 * m->__get(14))))->init(7,((((A_4 * m->__get(3)) + (A_5 * m->__get(7))) + (A_6 * m->__get(11))) + (A_7 * m->__get(15))))->init(8,((((A_8 * m->__get(0)) + (A_9 * m->__get(4))) + (A_10 * m->__get(8))) + (A_11 * m->__get(12))))->init(9,((((A_8 * m->__get(1)) + (A_9 * m->__get(5))) + (A_10 * m->__get(9))) + (A_11 * m->__get(13))))->init(10,((((A_8 * m->__get(2)) + (A_9 * m->__get(6))) + (A_10 * m->__get(10))) + (A_11 * m->__get(14))))->init(11,((((A_8 * m->__get(3)) + (A_9 * m->__get(7))) + (A_10 * m->__get(11))) + (A_11 * m->__get(15))))->init(12,((((A_12 * m->__get(0)) + (A_13 * m->__get(4))) + (A_14 * m->__get(8))) + (A_15 * m->__get(12))))->init(13,((((A_12 * m->__get(1)) + (A_13 * m->__get(5))) + (A_14 * m->__get(9))) + (A_15 * m->__get(13))))->init(14,((((A_12 * m->__get(2)) + (A_13 * m->__get(6))) + (A_14 * m->__get(10))) + (A_15 * m->__get(14))))->init(15,((((A_12 * m->__get(3)) + (A_13 * m->__get(7))) + (A_14 * m->__get(11))) + (A_15 * m->__get(15))));
            							}
            							else {
HXLINE( 138)								if (::StringTools_obj::startsWith(t,HX_("rotateZ",9f,49,1d,f1))) {
HXLINE( 139)									Float a = ((nums->__get(0) * ::Math_obj::PI) / ( (Float)(180) ));
HXDLIN( 139)									Float A_0 = ::Math_obj::cos(a);
HXDLIN( 139)									Float A_1 = -(::Math_obj::sin(a));
HXDLIN( 139)									Float A_2 = ( (Float)(0) );
HXDLIN( 139)									Float A_3 = ( (Float)(0) );
HXDLIN( 139)									Float A_4 = ::Math_obj::sin(a);
HXDLIN( 139)									Float A_5 = ::Math_obj::cos(a);
HXDLIN( 139)									Float A_6 = ( (Float)(0) );
HXDLIN( 139)									Float A_7 = ( (Float)(0) );
HXDLIN( 139)									Float A_8 = ( (Float)(0) );
HXDLIN( 139)									Float A_9 = ( (Float)(0) );
HXDLIN( 139)									Float A_10 = ( (Float)(1) );
HXDLIN( 139)									Float A_11 = ( (Float)(0) );
HXDLIN( 139)									Float A_12 = ( (Float)(0) );
HXDLIN( 139)									Float A_13 = ( (Float)(0) );
HXDLIN( 139)									Float A_14 = ( (Float)(0) );
HXDLIN( 139)									Float A_15 = ( (Float)(1) );
HXDLIN( 139)									m = ::Array_obj< Float >::__new(16)->init(0,((((A_0 * m->__get(0)) + (A_1 * m->__get(4))) + (A_2 * m->__get(8))) + (A_3 * m->__get(12))))->init(1,((((A_0 * m->__get(1)) + (A_1 * m->__get(5))) + (A_2 * m->__get(9))) + (A_3 * m->__get(13))))->init(2,((((A_0 * m->__get(2)) + (A_1 * m->__get(6))) + (A_2 * m->__get(10))) + (A_3 * m->__get(14))))->init(3,((((A_0 * m->__get(3)) + (A_1 * m->__get(7))) + (A_2 * m->__get(11))) + (A_3 * m->__get(15))))->init(4,((((A_4 * m->__get(0)) + (A_5 * m->__get(4))) + (A_6 * m->__get(8))) + (A_7 * m->__get(12))))->init(5,((((A_4 * m->__get(1)) + (A_5 * m->__get(5))) + (A_6 * m->__get(9))) + (A_7 * m->__get(13))))->init(6,((((A_4 * m->__get(2)) + (A_5 * m->__get(6))) + (A_6 * m->__get(10))) + (A_7 * m->__get(14))))->init(7,((((A_4 * m->__get(3)) + (A_5 * m->__get(7))) + (A_6 * m->__get(11))) + (A_7 * m->__get(15))))->init(8,((((A_8 * m->__get(0)) + (A_9 * m->__get(4))) + (A_10 * m->__get(8))) + (A_11 * m->__get(12))))->init(9,((((A_8 * m->__get(1)) + (A_9 * m->__get(5))) + (A_10 * m->__get(9))) + (A_11 * m->__get(13))))->init(10,((((A_8 * m->__get(2)) + (A_9 * m->__get(6))) + (A_10 * m->__get(10))) + (A_11 * m->__get(14))))->init(11,((((A_8 * m->__get(3)) + (A_9 * m->__get(7))) + (A_10 * m->__get(11))) + (A_11 * m->__get(15))))->init(12,((((A_12 * m->__get(0)) + (A_13 * m->__get(4))) + (A_14 * m->__get(8))) + (A_15 * m->__get(12))))->init(13,((((A_12 * m->__get(1)) + (A_13 * m->__get(5))) + (A_14 * m->__get(9))) + (A_15 * m->__get(13))))->init(14,((((A_12 * m->__get(2)) + (A_13 * m->__get(6))) + (A_14 * m->__get(10))) + (A_15 * m->__get(14))))->init(15,((((A_12 * m->__get(3)) + (A_13 * m->__get(7))) + (A_14 * m->__get(11))) + (A_15 * m->__get(15))));
            								}
            								else {
HXLINE( 140)									if (::StringTools_obj::startsWith(t,HX_("scale",8a,ce,ce,78))) {
HXLINE( 141)										Float A_0 = nums->__get(0);
HXDLIN( 141)										Float A_1 = ( (Float)(0) );
HXDLIN( 141)										Float A_2 = ( (Float)(0) );
HXDLIN( 141)										Float A_3 = ( (Float)(0) );
HXDLIN( 141)										Float A_4 = ( (Float)(0) );
HXDLIN( 141)										Float A_5 = nums->__get(1);
HXDLIN( 141)										Float A_6 = ( (Float)(0) );
HXDLIN( 141)										Float A_7 = ( (Float)(0) );
HXDLIN( 141)										Float A_8 = ( (Float)(0) );
HXDLIN( 141)										Float A_9 = ( (Float)(0) );
HXDLIN( 141)										Float A_10 = nums->__get(2);
HXDLIN( 141)										Float A_11 = ( (Float)(0) );
HXDLIN( 141)										Float A_12 = ( (Float)(0) );
HXDLIN( 141)										Float A_13 = ( (Float)(0) );
HXDLIN( 141)										Float A_14 = ( (Float)(0) );
HXDLIN( 141)										Float A_15 = ( (Float)(1) );
HXDLIN( 141)										m = ::Array_obj< Float >::__new(16)->init(0,((((A_0 * m->__get(0)) + (A_1 * m->__get(4))) + (A_2 * m->__get(8))) + (A_3 * m->__get(12))))->init(1,((((A_0 * m->__get(1)) + (A_1 * m->__get(5))) + (A_2 * m->__get(9))) + (A_3 * m->__get(13))))->init(2,((((A_0 * m->__get(2)) + (A_1 * m->__get(6))) + (A_2 * m->__get(10))) + (A_3 * m->__get(14))))->init(3,((((A_0 * m->__get(3)) + (A_1 * m->__get(7))) + (A_2 * m->__get(11))) + (A_3 * m->__get(15))))->init(4,((((A_4 * m->__get(0)) + (A_5 * m->__get(4))) + (A_6 * m->__get(8))) + (A_7 * m->__get(12))))->init(5,((((A_4 * m->__get(1)) + (A_5 * m->__get(5))) + (A_6 * m->__get(9))) + (A_7 * m->__get(13))))->init(6,((((A_4 * m->__get(2)) + (A_5 * m->__get(6))) + (A_6 * m->__get(10))) + (A_7 * m->__get(14))))->init(7,((((A_4 * m->__get(3)) + (A_5 * m->__get(7))) + (A_6 * m->__get(11))) + (A_7 * m->__get(15))))->init(8,((((A_8 * m->__get(0)) + (A_9 * m->__get(4))) + (A_10 * m->__get(8))) + (A_11 * m->__get(12))))->init(9,((((A_8 * m->__get(1)) + (A_9 * m->__get(5))) + (A_10 * m->__get(9))) + (A_11 * m->__get(13))))->init(10,((((A_8 * m->__get(2)) + (A_9 * m->__get(6))) + (A_10 * m->__get(10))) + (A_11 * m->__get(14))))->init(11,((((A_8 * m->__get(3)) + (A_9 * m->__get(7))) + (A_10 * m->__get(11))) + (A_11 * m->__get(15))))->init(12,((((A_12 * m->__get(0)) + (A_13 * m->__get(4))) + (A_14 * m->__get(8))) + (A_15 * m->__get(12))))->init(13,((((A_12 * m->__get(1)) + (A_13 * m->__get(5))) + (A_14 * m->__get(9))) + (A_15 * m->__get(13))))->init(14,((((A_12 * m->__get(2)) + (A_13 * m->__get(6))) + (A_14 * m->__get(10))) + (A_15 * m->__get(14))))->init(15,((((A_12 * m->__get(3)) + (A_13 * m->__get(7))) + (A_14 * m->__get(11))) + (A_15 * m->__get(15))));
            									}
            									else {
HXLINE( 142)										if (::StringTools_obj::startsWith(t,HX_("translate",4e,d7,7f,49))) {
HXLINE( 143)											Float A_0 = ( (Float)(1) );
HXDLIN( 143)											Float A_1 = ( (Float)(0) );
HXDLIN( 143)											Float A_2 = ( (Float)(0) );
HXDLIN( 143)											Float A_3 = nums->__get(0);
HXDLIN( 143)											Float A_4 = ( (Float)(0) );
HXDLIN( 143)											Float A_5 = ( (Float)(1) );
HXDLIN( 143)											Float A_6 = ( (Float)(0) );
HXDLIN( 143)											Float A_7 = nums->__get(1);
HXDLIN( 143)											Float A_8 = ( (Float)(0) );
HXDLIN( 143)											Float A_9 = ( (Float)(0) );
HXDLIN( 143)											Float A_10 = ( (Float)(1) );
HXDLIN( 143)											Float A_11 = nums->__get(2);
HXDLIN( 143)											Float A_12 = ( (Float)(0) );
HXDLIN( 143)											Float A_13 = ( (Float)(0) );
HXDLIN( 143)											Float A_14 = ( (Float)(0) );
HXDLIN( 143)											Float A_15 = ( (Float)(1) );
HXDLIN( 143)											m = ::Array_obj< Float >::__new(16)->init(0,((((A_0 * m->__get(0)) + (A_1 * m->__get(4))) + (A_2 * m->__get(8))) + (A_3 * m->__get(12))))->init(1,((((A_0 * m->__get(1)) + (A_1 * m->__get(5))) + (A_2 * m->__get(9))) + (A_3 * m->__get(13))))->init(2,((((A_0 * m->__get(2)) + (A_1 * m->__get(6))) + (A_2 * m->__get(10))) + (A_3 * m->__get(14))))->init(3,((((A_0 * m->__get(3)) + (A_1 * m->__get(7))) + (A_2 * m->__get(11))) + (A_3 * m->__get(15))))->init(4,((((A_4 * m->__get(0)) + (A_5 * m->__get(4))) + (A_6 * m->__get(8))) + (A_7 * m->__get(12))))->init(5,((((A_4 * m->__get(1)) + (A_5 * m->__get(5))) + (A_6 * m->__get(9))) + (A_7 * m->__get(13))))->init(6,((((A_4 * m->__get(2)) + (A_5 * m->__get(6))) + (A_6 * m->__get(10))) + (A_7 * m->__get(14))))->init(7,((((A_4 * m->__get(3)) + (A_5 * m->__get(7))) + (A_6 * m->__get(11))) + (A_7 * m->__get(15))))->init(8,((((A_8 * m->__get(0)) + (A_9 * m->__get(4))) + (A_10 * m->__get(8))) + (A_11 * m->__get(12))))->init(9,((((A_8 * m->__get(1)) + (A_9 * m->__get(5))) + (A_10 * m->__get(9))) + (A_11 * m->__get(13))))->init(10,((((A_8 * m->__get(2)) + (A_9 * m->__get(6))) + (A_10 * m->__get(10))) + (A_11 * m->__get(14))))->init(11,((((A_8 * m->__get(3)) + (A_9 * m->__get(7))) + (A_10 * m->__get(11))) + (A_11 * m->__get(15))))->init(12,((((A_12 * m->__get(0)) + (A_13 * m->__get(4))) + (A_14 * m->__get(8))) + (A_15 * m->__get(12))))->init(13,((((A_12 * m->__get(1)) + (A_13 * m->__get(5))) + (A_14 * m->__get(9))) + (A_15 * m->__get(13))))->init(14,((((A_12 * m->__get(2)) + (A_13 * m->__get(6))) + (A_14 * m->__get(10))) + (A_15 * m->__get(14))))->init(15,((((A_12 * m->__get(3)) + (A_13 * m->__get(7))) + (A_14 * m->__get(11))) + (A_15 * m->__get(15))));
            										}
            										else {
HXLINE( 144)											if (::StringTools_obj::startsWith(t,HX_("matrix",41,36,c8,bb))) {
HXLINE( 145)												m = ::Array_obj< Float >::__new(16)->init(0,((((nums->__get(0) * m->__get(0)) + (nums->__get(1) * m->__get(4))) + (nums->__get(2) * m->__get(8))) + (nums->__get(3) * m->__get(12))))->init(1,((((nums->__get(0) * m->__get(1)) + (nums->__get(1) * m->__get(5))) + (nums->__get(2) * m->__get(9))) + (nums->__get(3) * m->__get(13))))->init(2,((((nums->__get(0) * m->__get(2)) + (nums->__get(1) * m->__get(6))) + (nums->__get(2) * m->__get(10))) + (nums->__get(3) * m->__get(14))))->init(3,((((nums->__get(0) * m->__get(3)) + (nums->__get(1) * m->__get(7))) + (nums->__get(2) * m->__get(11))) + (nums->__get(3) * m->__get(15))))->init(4,((((nums->__get(4) * m->__get(0)) + (nums->__get(5) * m->__get(4))) + (nums->__get(6) * m->__get(8))) + (nums->__get(7) * m->__get(12))))->init(5,((((nums->__get(4) * m->__get(1)) + (nums->__get(5) * m->__get(5))) + (nums->__get(6) * m->__get(9))) + (nums->__get(7) * m->__get(13))))->init(6,((((nums->__get(4) * m->__get(2)) + (nums->__get(5) * m->__get(6))) + (nums->__get(6) * m->__get(10))) + (nums->__get(7) * m->__get(14))))->init(7,((((nums->__get(4) * m->__get(3)) + (nums->__get(5) * m->__get(7))) + (nums->__get(6) * m->__get(11))) + (nums->__get(7) * m->__get(15))))->init(8,((((nums->__get(8) * m->__get(0)) + (nums->__get(9) * m->__get(4))) + (nums->__get(10) * m->__get(8))) + (nums->__get(11) * m->__get(12))))->init(9,((((nums->__get(8) * m->__get(1)) + (nums->__get(9) * m->__get(5))) + (nums->__get(10) * m->__get(9))) + (nums->__get(11) * m->__get(13))))->init(10,((((nums->__get(8) * m->__get(2)) + (nums->__get(9) * m->__get(6))) + (nums->__get(10) * m->__get(10))) + (nums->__get(11) * m->__get(14))))->init(11,((((nums->__get(8) * m->__get(3)) + (nums->__get(9) * m->__get(7))) + (nums->__get(10) * m->__get(11))) + (nums->__get(11) * m->__get(15))))->init(12,((((nums->__get(12) * m->__get(0)) + (nums->__get(13) * m->__get(4))) + (nums->__get(14) * m->__get(8))) + (nums->__get(15) * m->__get(12))))->init(13,((((nums->__get(12) * m->__get(1)) + (nums->__get(13) * m->__get(5))) + (nums->__get(14) * m->__get(9))) + (nums->__get(15) * m->__get(13))))->init(14,((((nums->__get(12) * m->__get(2)) + (nums->__get(13) * m->__get(6))) + (nums->__get(14) * m->__get(10))) + (nums->__get(15) * m->__get(14))))->init(15,((((nums->__get(12) * m->__get(3)) + (nums->__get(13) * m->__get(7))) + (nums->__get(14) * m->__get(11))) + (nums->__get(15) * m->__get(15))));
            											}
            											else {
HXLINE( 146)												if (::StringTools_obj::startsWith(t,HX_("focal",65,4e,89,04))) {
HXLINE( 147)													render1->setFocal(nums->__get(0));
            												}
            												else {
HXLINE( 148)													if (::StringTools_obj::startsWith(t,HX_("auto",6f,df,76,40))) {
HXLINE( 149)														render1->transform(m);
HXLINE( 150)														m = ::Array_obj< Float >::fromData( _hx_array_data_00333580_9,16);
HXLINE( 151)														if ((nums->length == 0)) {
HXLINE( 152)															render1->autoPlace(null(),null());
            														}
            														else {
HXLINE( 154)															render1->autoPlace(nums->__get(0),nums->__get(1));
            														}
            													}
            													else {
HXLINE( 157)														::haxe::Log_obj::trace((HX_("[warn] invalid transform expression: ",ff,93,6b,3d) + t),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),157,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            													}
            												}
            											}
            										}
            									}
            								}
            							}
            						}
            					}
            				}
HXLINE( 160)				render1->transform(m);
HXLINE( 163)				bool _hx_tmp;
HXDLIN( 163)				if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) != HX_("edges",96,6d,e0,69))) {
HXLINE( 163)					_hx_tmp = (kwargs->get_string(HX_("--mode",c3,2a,78,2f)) != HX_("vertices",f9,bf,15,6a));
            				}
            				else {
HXLINE( 163)					_hx_tmp = false;
            				}
HXDLIN( 163)				if (_hx_tmp) {
HXLINE( 164)					if ((kwargs->get_string(HX_("--cull",12,9d,e0,28)) == HX_("true",4e,a7,03,4d))) {
HXLINE( 165)						render1->apparentRidges(::Std_obj::parseFloat(kwargs->get_string(HX_("--thresh",dc,79,dd,eb))),null());
            					}
            					else {
HXLINE( 166)						if ((kwargs->get_string(HX_("--cull",12,9d,e0,28)) == HX_("false",a3,35,4f,fb))) {
HXLINE( 167)							render1->apparentRidges(::Std_obj::parseFloat(kwargs->get_string(HX_("--thresh",dc,79,dd,eb))),-1);
            						}
            						else {
HXLINE( 169)							Float _hx_tmp = ::Std_obj::parseFloat(kwargs->get_string(HX_("--thresh",dc,79,dd,eb)));
HXDLIN( 169)							render1->apparentRidges(_hx_tmp,::Std_obj::parseFloat(kwargs->get_string(HX_("--cull",12,9d,e0,28))));
            						}
            					}
            				}
HXLINE( 173)				if ((verb > 0)) {
HXLINE( 173)					::haxe::Log_obj::trace(HX_("[info] generating output...",4b,3b,d9,dd),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),173,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            				}
HXLINE( 175)				::String out = HX_("",00,00,00,00);
HXLINE( 177)				if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("edges",96,6d,e0,69))) {
HXLINE( 178)					render1->edges();
HXLINE( 179)					out = ::apparentridges::SVGWriter_obj::lines(render1,null());
            				}
            				else {
HXLINE( 180)					if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("vertices",f9,bf,15,6a))) {
HXLINE( 181)						render1->vertices();
HXLINE( 182)						out = ::apparentridges::SVGWriter_obj::lines(render1,null());
            					}
            					else {
HXLINE( 183)						if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("lines",ff,dd,01,75))) {
HXLINE( 184)							out = ::apparentridges::SVGWriter_obj::lines(render1,null());
            						}
            						else {
HXLINE( 185)							if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("polylines",33,0c,bc,77))) {
HXLINE( 186)								render1->buildPolylines(null());
HXLINE( 187)								out = ::apparentridges::SVGWriter_obj::polylines(render1,null());
            							}
            							else {
HXLINE( 188)								if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("gradients",83,78,13,cd))) {
HXLINE( 189)									render1->buildPolylines(null());
HXLINE( 190)									out = ::apparentridges::SVGWriter_obj::gradients(render1,null());
            								}
            								else {
HXLINE( 191)									if ((kwargs->get_string(HX_("--mode",c3,2a,78,2f)) == HX_("pixels",2d,ef,a9,8c))) {
HXLINE( 195)										if ((verb > 0)) {
HXLINE( 195)											::haxe::Log_obj::trace(HX_("[info] making depth map...",b8,56,ee,eb),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),195,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 192)										::Array< Float > d = ::apparentridges::PixelMap_obj::depth(render1,true);
HXLINE( 197)										if ((verb > 0)) {
HXLINE( 197)											::haxe::Log_obj::trace(HX_("[info] making ppm string...",c1,b5,65,d8),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),197,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 193)										::String o = ::apparentridges::PixelMap_obj::toPPMString(d,render1->width,render1->height,( (Float)(0) ),( (Float)(1) ));
HXLINE( 194)										::String p = (inp + HX_(".depth.ppm",94,e7,c8,17));
HXLINE( 200)										::sys::io::File_obj::saveContent(p,o);
HXLINE( 201)										if ((verb > 0)) {
HXLINE( 201)											::haxe::Log_obj::trace((HX_("[info] wrote ",c7,5e,ba,b2) + p),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),201,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 203)										if ((verb > 0)) {
HXLINE( 203)											::haxe::Log_obj::trace(HX_("[info] making normal map...",22,49,a1,98),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),203,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 204)										d = ::apparentridges::PixelMap_obj::normal(render1);
HXLINE( 205)										o = ::apparentridges::PixelMap_obj::toPPMString(d,render1->width,render1->height,( (Float)(-1) ),( (Float)(1) ));
HXLINE( 206)										p = (inp + HX_(".normal.ppm",54,86,5f,82));
HXLINE( 207)										::sys::io::File_obj::saveContent(p,o);
HXLINE( 208)										if ((verb > 0)) {
HXLINE( 208)											::haxe::Log_obj::trace((HX_("[info] wrote ",c7,5e,ba,b2) + p),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),208,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 210)										if ((verb > 0)) {
HXLINE( 210)											::haxe::Log_obj::trace(HX_("[info] making curvature map...",1c,52,cb,19),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),210,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 211)										d = ::apparentridges::PixelMap_obj::curvature(render1);
HXLINE( 212)										o = ::apparentridges::PixelMap_obj::toPPMString(d,render1->width,render1->height,((Float)-0.5),((Float)0.5));
HXLINE( 213)										p = (inp + HX_(".curv.ppm",23,0b,e7,ab));
HXLINE( 214)										::sys::io::File_obj::saveContent(p,o);
HXLINE( 215)										if ((verb > 0)) {
HXLINE( 215)											::haxe::Log_obj::trace((HX_("[info] wrote ",c7,5e,ba,b2) + p),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),215,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            										}
HXLINE( 217)										continue;
            									}
            									else {
HXLINE( 220)										::haxe::Log_obj::trace(HX_("[error] invalid mode",1c,c1,54,d3),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),220,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
HXLINE( 221)										return 1;
            									}
            								}
            							}
            						}
            					}
            				}
HXLINE( 223)				if ((verb > 0)) {
HXLINE( 223)					::haxe::Log_obj::trace(HX_("[info] writing file...",82,4c,36,f4),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),223,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            				}
HXLINE( 225)				::String path = kwargs->get_string(HX_("--output",61,fe,a4,69));
HXLINE( 226)				if ((path == HX_("stdout",cb,bf,f3,07))) {
HXLINE( 227)					::Sys_obj::println(out);
            				}
            				else {
HXLINE( 229)					if ((args->length > 1)) {
HXLINE( 230)						path = (inp + HX_(".svg",f6,7a,bf,1e));
            					}
HXLINE( 232)					::sys::io::File_obj::saveContent(path,out);
HXLINE( 233)					if ((verb > 0)) {
HXLINE( 233)						::haxe::Log_obj::trace((HX_("[info] wrote ",c7,5e,ba,b2) + path),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),233,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            					}
            				}
            			}
            		}
HXLINE( 237)		if ((verb > 0)) {
HXLINE( 237)			::haxe::Log_obj::trace(HX_("[info] all files processed.",74,0e,5a,c3),::hx::SourceInfo(HX_("Cli.hx",3e,7d,ba,3e),237,HX_("Cli",80,35,33,00),HX_("main",39,38,56,48)));
            		}
HXLINE( 238)		return 0;
            	}


STATIC_HX_DEFINE_DYNAMIC_FUNC0(Cli_obj,main,return )


Cli_obj::Cli_obj()
{
}

bool Cli_obj::__GetStatic(const ::String &inName, Dynamic &outValue, ::hx::PropertyAccess inCallProp)
{
	switch(inName.length) {
	case 4:
		if (HX_FIELD_EQ(inName,"main") ) { outValue = main_dyn(); return true; }
	}
	return false;
}

#ifdef HXCPP_SCRIPTABLE
static ::hx::StorageInfo *Cli_obj_sMemberStorageInfo = 0;
static ::hx::StaticInfo *Cli_obj_sStaticStorageInfo = 0;
#endif

::hx::Class Cli_obj::__mClass;

static ::String Cli_obj_sStaticFields[] = {
	HX_("main",39,38,56,48),
	::String(null())
};

void Cli_obj::__register()
{
	Cli_obj _hx_dummy;
	Cli_obj::_hx_vtable = *(void **)&_hx_dummy;
	::hx::Static(__mClass) = new ::hx::Class_obj();
	__mClass->mName = HX_("Cli",80,35,33,00);
	__mClass->mSuper = &super::__SGetClass();
	__mClass->mConstructEmpty = &__CreateEmpty;
	__mClass->mConstructArgs = &__Create;
	__mClass->mGetStaticField = &Cli_obj::__GetStatic;
	__mClass->mSetStaticField = &::hx::Class_obj::SetNoStaticField;
	__mClass->mStatics = ::hx::Class_obj::dupFunctions(Cli_obj_sStaticFields);
	__mClass->mMembers = ::hx::Class_obj::dupFunctions(0 /* sMemberFields */);
	__mClass->mCanCast = ::hx::TCanCast< Cli_obj >;
#ifdef HXCPP_SCRIPTABLE
	__mClass->mMemberStorageInfo = Cli_obj_sMemberStorageInfo;
#endif
#ifdef HXCPP_SCRIPTABLE
	__mClass->mStaticStorageInfo = Cli_obj_sStaticStorageInfo;
#endif
	::hx::_hx_RegisterClass(__mClass->mName, __mClass);
}

