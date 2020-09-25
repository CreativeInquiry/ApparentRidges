// Generated by Haxe 4.1.3
#include <hxcpp.h>

#ifndef INCLUDED_95f339a1d026d52c
#define INCLUDED_95f339a1d026d52c
#include "hxMath.h"
#endif
#ifndef INCLUDED_Std
#include <Std.h>
#endif
#ifndef INCLUDED_StringBuf
#include <StringBuf.h>
#endif
#ifndef INCLUDED_apparentridges_Line
#include <apparentridges/Line.h>
#endif
#ifndef INCLUDED_apparentridges_Render
#include <apparentridges/Render.h>
#endif
#ifndef INCLUDED_apparentridges_SVGWriter
#include <apparentridges/SVGWriter.h>
#endif

HX_LOCAL_STACK_FRAME(_hx_pos_7941ea06c644aa35_1881_rd,"apparentridges.SVGWriter","rd",0x7e7f1636,"apparentridges.SVGWriter.rd","apparentridges/ApparentRidges.hx",1881,0xfeacc84f)
HX_LOCAL_STACK_FRAME(_hx_pos_7941ea06c644aa35_1884_lines,"apparentridges.SVGWriter","lines",0x987e95fb,"apparentridges.SVGWriter.lines","apparentridges/ApparentRidges.hx",1884,0xfeacc84f)
static const ::String _hx_array_data_04a6b36a_3[] = {
	HX_("</svg>\n",1d,0c,30,0f),
};
HX_LOCAL_STACK_FRAME(_hx_pos_7941ea06c644aa35_1903_polylines,"apparentridges.SVGWriter","polylines",0x765a722f,"apparentridges.SVGWriter.polylines","apparentridges/ApparentRidges.hx",1903,0xfeacc84f)
static const ::String _hx_array_data_04a6b36a_7[] = {
	HX_("  <polyline points=\"",4c,5a,8f,f0),
};
static const ::String _hx_array_data_04a6b36a_8[] = {
	HX_("</svg>\n",1d,0c,30,0f),
};
HX_LOCAL_STACK_FRAME(_hx_pos_7941ea06c644aa35_1925_gradients,"apparentridges.SVGWriter","gradients",0xcbb1de7f,"apparentridges.SVGWriter.gradients","apparentridges/ApparentRidges.hx",1925,0xfeacc84f)
static const ::String _hx_array_data_04a6b36a_15[] = {
	HX_("  <g fill=\"none\" stroke-width=\"1\">\n",66,1e,73,ea),
};
static const ::String _hx_array_data_04a6b36a_16[] = {
	HX_("  </g>\n",a0,97,84,a9),
};
static const ::String _hx_array_data_04a6b36a_17[] = {
	HX_("</svg>\n",1d,0c,30,0f),
};
namespace apparentridges{

void SVGWriter_obj::__construct() { }

Dynamic SVGWriter_obj::__CreateEmpty() { return new SVGWriter_obj; }

void *SVGWriter_obj::_hx_vtable = 0;

Dynamic SVGWriter_obj::__Create(::hx::DynamicArray inArgs)
{
	::hx::ObjectPtr< SVGWriter_obj > _hx_result = new SVGWriter_obj();
	_hx_result->__construct();
	return _hx_result;
}

bool SVGWriter_obj::_hx_isInstanceOf(int inClassId) {
	return inClassId==(int)0x00000001 || inClassId==(int)0x48ccdc52;
}

Float SVGWriter_obj::rd(Float x){
            	HX_STACKFRAME(&_hx_pos_7941ea06c644aa35_1881_rd)
HXDLIN(1881)		return (( (Float)(::Math_obj::round((x * ( (Float)(100) )))) ) / ( (Float)(100) ));
            	}


STATIC_HX_DEFINE_DYNAMIC_FUNC1(SVGWriter_obj,rd,return )

::String SVGWriter_obj::lines( ::apparentridges::Render render,::hx::Null< bool >  __o_useOpacity){
            		bool useOpacity = __o_useOpacity.Default(true);
            	HX_GC_STACKFRAME(&_hx_pos_7941ea06c644aa35_1884_lines)
HXLINE(1885)		int w = render->width;
HXLINE(1886)		int h = render->height;
HXLINE(1887)		 ::StringBuf out =  ::StringBuf_obj::__alloc( HX_CTX );
HXLINE(1888)		{
HXLINE(1888)			::String x = ((((HX_("<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"",d2,21,50,12) + w) + HX_("\" height=\"",aa,50,5c,1b)) + h) + HX_("\">\n",ae,02,1a,00));
HXDLIN(1888)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1888)				out->flush();
            			}
HXDLIN(1888)			if (::hx::IsNull( out->b )) {
HXLINE(1888)				out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x));
            			}
            			else {
HXLINE(1888)				::Array< ::String > out1 = out->b;
HXDLIN(1888)				out1->push(::Std_obj::string(x));
            			}
            		}
HXLINE(1890)		{
HXLINE(1890)			int _g = 0;
HXDLIN(1890)			int _g1 = render->lines->length;
HXDLIN(1890)			while((_g < _g1)){
HXLINE(1890)				_g = (_g + 1);
HXDLIN(1890)				int i = (_g - 1);
HXLINE(1891)				Float x1 = (( (Float)(::Math_obj::round((render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->x1 * ( (Float)(100) )))) ) / ( (Float)(100) ));
HXLINE(1892)				Float y1 = (( (Float)(::Math_obj::round((render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->y1 * ( (Float)(100) )))) ) / ( (Float)(100) ));
HXLINE(1893)				Float x2 = (( (Float)(::Math_obj::round((render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->x2 * ( (Float)(100) )))) ) / ( (Float)(100) ));
HXLINE(1894)				Float y2 = (( (Float)(::Math_obj::round((render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->y2 * ( (Float)(100) )))) ) / ( (Float)(100) ));
HXLINE(1895)				Float o;
HXDLIN(1895)				if (useOpacity) {
HXLINE(1895)					o = ((render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->opacity1 + render->lines->__get(i).StaticCast<  ::apparentridges::Line >()->opacity2) / ( (Float)(2) ));
            				}
            				else {
HXLINE(1895)					o = ( (Float)(1) );
            				}
HXLINE(1896)				int oi = ::Std_obj::_hx_int((( (Float)(255) ) - (o * ( (Float)(255) ))));
HXLINE(1897)				{
HXLINE(1897)					::String x = ((((((((((((((HX_("  <line x1=\"",2e,56,46,d5) + x1) + HX_("\" y1=\"",9b,09,15,55)) + y1) + HX_("\" x2=\"",3d,95,6c,54)) + x2) + HX_("\" y2=\"",dc,cb,15,55)) + y2) + HX_("\" fill=\"none\" stroke=\"rgb(",b4,15,9e,c9)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(")\" stroke-width=\"1\" stroke-linecap=\"round\"/>\n",27,4b,34,0c));
HXDLIN(1897)					if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1897)						out->flush();
            					}
HXDLIN(1897)					if (::hx::IsNull( out->b )) {
HXLINE(1897)						out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x));
            					}
            					else {
HXLINE(1897)						::Array< ::String > out1 = out->b;
HXDLIN(1897)						out1->push(::Std_obj::string(x));
            					}
            				}
            			}
            		}
HXLINE(1899)		{
HXLINE(1899)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1899)				out->flush();
            			}
HXDLIN(1899)			if (::hx::IsNull( out->b )) {
HXLINE(1899)				out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_3,1);
            			}
            			else {
HXLINE(1899)				out->b->push(HX_("</svg>\n",1d,0c,30,0f));
            			}
            		}
HXLINE(1900)		return out->toString();
            	}


STATIC_HX_DEFINE_DYNAMIC_FUNC2(SVGWriter_obj,lines,return )

::String SVGWriter_obj::polylines( ::apparentridges::Render render,::hx::Null< bool >  __o_colorful){
            		bool colorful = __o_colorful.Default(false);
            	HX_GC_STACKFRAME(&_hx_pos_7941ea06c644aa35_1903_polylines)
HXLINE(1904)		int w = render->width;
HXLINE(1905)		int h = render->height;
HXLINE(1906)		 ::StringBuf out =  ::StringBuf_obj::__alloc( HX_CTX );
HXLINE(1907)		{
HXLINE(1907)			::String x = ((((HX_("<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"",d2,21,50,12) + w) + HX_("\" height=\"",aa,50,5c,1b)) + h) + HX_("\">\n",ae,02,1a,00));
HXDLIN(1907)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1907)				out->flush();
            			}
HXDLIN(1907)			if (::hx::IsNull( out->b )) {
HXLINE(1907)				out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x));
            			}
            			else {
HXLINE(1907)				::Array< ::String > out1 = out->b;
HXDLIN(1907)				out1->push(::Std_obj::string(x));
            			}
            		}
HXLINE(1909)		::String color = HX_("black",bf,d5,f1,b4);
HXLINE(1910)		{
HXLINE(1910)			int _g = 0;
HXDLIN(1910)			int _g1 = render->polylines->length;
HXDLIN(1910)			while((_g < _g1)){
HXLINE(1910)				_g = (_g + 1);
HXDLIN(1910)				int i = (_g - 1);
HXLINE(1911)				{
HXLINE(1911)					if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1911)						out->flush();
            					}
HXDLIN(1911)					if (::hx::IsNull( out->b )) {
HXLINE(1911)						out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_7,1);
            					}
            					else {
HXLINE(1911)						out->b->push(HX_("  <polyline points=\"",4c,5a,8f,f0));
            					}
            				}
HXLINE(1912)				{
HXLINE(1912)					int _g1 = 0;
HXDLIN(1912)					int _g2 = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->length;
HXDLIN(1912)					while((_g1 < _g2)){
HXLINE(1912)						_g1 = (_g1 + 1);
HXDLIN(1912)						int j = (_g1 - 1);
HXLINE(1913)						::Array< Float > p = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->__get(j).StaticCast< ::Array< Float > >();
HXLINE(1914)						{
HXLINE(1914)							::String x = ((HX_("",00,00,00,00) + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(p,0)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_(",",2c,00,00,00));
HXDLIN(1914)							::String x1 = ((x + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(p,1)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_(" ",20,00,00,00));
HXDLIN(1914)							if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1914)								out->flush();
            							}
HXDLIN(1914)							if (::hx::IsNull( out->b )) {
HXLINE(1914)								out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x1));
            							}
            							else {
HXLINE(1914)								::Array< ::String > out1 = out->b;
HXDLIN(1914)								out1->push(::Std_obj::string(x1));
            							}
            						}
            					}
            				}
HXLINE(1916)				if (colorful) {
HXLINE(1917)					::String color1 = ((HX_("rgb(",7b,d0,a8,4b) + ::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(128) )))) + HX_(",",2c,00,00,00));
HXDLIN(1917)					::String color2 = ((color1 + ::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(128) )))) + HX_(",",2c,00,00,00));
HXDLIN(1917)					color = ((color2 + ::Std_obj::_hx_int((::Math_obj::random() * ( (Float)(128) )))) + HX_(")",29,00,00,00));
            				}
HXLINE(1919)				{
HXLINE(1919)					::String x = ((HX_("\" fill=\"none\" stroke=\"",b9,f4,c9,9a) + color) + HX_("\" stroke-width=\"1\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n",b2,63,cb,91));
HXDLIN(1919)					if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1919)						out->flush();
            					}
HXDLIN(1919)					if (::hx::IsNull( out->b )) {
HXLINE(1919)						out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x));
            					}
            					else {
HXLINE(1919)						::Array< ::String > out1 = out->b;
HXDLIN(1919)						out1->push(::Std_obj::string(x));
            					}
            				}
            			}
            		}
HXLINE(1921)		{
HXLINE(1921)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1921)				out->flush();
            			}
HXDLIN(1921)			if (::hx::IsNull( out->b )) {
HXLINE(1921)				out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_8,1);
            			}
            			else {
HXLINE(1921)				out->b->push(HX_("</svg>\n",1d,0c,30,0f));
            			}
            		}
HXLINE(1922)		return out->toString();
            	}


STATIC_HX_DEFINE_DYNAMIC_FUNC2(SVGWriter_obj,polylines,return )

::String SVGWriter_obj::gradients( ::apparentridges::Render render,::hx::Null< Float >  __o_acc){
            		Float acc = __o_acc.Default(1);
            	HX_GC_STACKFRAME(&_hx_pos_7941ea06c644aa35_1925_gradients)
HXLINE(1926)		int w = render->width;
HXLINE(1927)		int h = render->height;
HXLINE(1928)		 ::StringBuf out =  ::StringBuf_obj::__alloc( HX_CTX );
HXLINE(1929)		{
HXLINE(1929)			::String x = ((((HX_("<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"",d2,21,50,12) + w) + HX_("\" height=\"",aa,50,5c,1b)) + h) + HX_("\">\n",ae,02,1a,00));
HXDLIN(1929)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1929)				out->flush();
            			}
HXDLIN(1929)			if (::hx::IsNull( out->b )) {
HXLINE(1929)				out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x));
            			}
            			else {
HXLINE(1929)				::Array< ::String > out1 = out->b;
HXDLIN(1929)				out1->push(::Std_obj::string(x));
            			}
            		}
HXLINE(1930)		int gid = 0;
HXLINE(1931)		{
HXLINE(1931)			int _g = 0;
HXDLIN(1931)			int _g1 = render->polylines->length;
HXDLIN(1931)			while((_g < _g1)){
HXLINE(1931)				_g = (_g + 1);
HXDLIN(1931)				int i = (_g - 1);
HXLINE(1932)				{
HXLINE(1932)					int _g1 = 0;
HXDLIN(1932)					int _g2 = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->length;
HXDLIN(1932)					while((_g1 < _g2)){
HXLINE(1932)						_g1 = (_g1 + 1);
HXDLIN(1932)						int j = (_g1 - 1);
HXLINE(1933)						::Array< Float > p = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->__get(j).StaticCast< ::Array< Float > >();
HXLINE(1934)						int oi = ::Std_obj::_hx_int((( (Float)(255) ) - (( (Float)(_hx_array_unsafe_get(p,2)) ) * ( (Float)(255) ))));
HXLINE(1935)						{
HXLINE(1935)							::String x = ((HX_("  <circle cx=\"",8e,d4,e0,d4) + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(p,0)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" cy=\"",39,fa,c0,46));
HXDLIN(1935)							::String x1 = ((((((((x + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(p,1)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" r=\"0.5\" stroke=\"none\" fill=\"rgb(",12,8e,3c,4d)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(")\"/>\n",c2,e9,0a,b2));
HXDLIN(1935)							if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1935)								out->flush();
            							}
HXDLIN(1935)							if (::hx::IsNull( out->b )) {
HXLINE(1935)								out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x1));
            							}
            							else {
HXLINE(1935)								::Array< ::String > out1 = out->b;
HXDLIN(1935)								out1->push(::Std_obj::string(x1));
            							}
            						}
            					}
            				}
            			}
            		}
HXLINE(1938)		{
HXLINE(1938)			int _g2 = 0;
HXDLIN(1938)			int _g3 = render->polylines->length;
HXDLIN(1938)			while((_g2 < _g3)){
HXLINE(1938)				_g2 = (_g2 + 1);
HXDLIN(1938)				int i = (_g2 - 1);
HXLINE(1939)				{
HXLINE(1939)					if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1939)						out->flush();
            					}
HXDLIN(1939)					if (::hx::IsNull( out->b )) {
HXLINE(1939)						out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_15,1);
            					}
            					else {
HXLINE(1939)						out->b->push(HX_("  <g fill=\"none\" stroke-width=\"1\">\n",66,1e,73,ea));
            					}
            				}
HXLINE(1940)				{
HXLINE(1940)					int _g = 0;
HXDLIN(1940)					int _g1 = (render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->length - 1);
HXDLIN(1940)					while((_g < _g1)){
HXLINE(1940)						_g = (_g + 1);
HXDLIN(1940)						int j = (_g - 1);
HXLINE(1941)						::Array< Float > p = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->__get(j).StaticCast< ::Array< Float > >();
HXLINE(1942)						::Array< Float > q = render->polylines->__get(i).StaticCast< ::Array< ::Dynamic> >()->__get((j + 1)).StaticCast< ::Array< Float > >();
HXLINE(1943)						Float d = ( (Float)(_hx_array_unsafe_get(p,2)) );
HXDLIN(1943)						int d1 = (::Std_obj::_hx_int(((::Math_obj::abs((d - ( (Float)(_hx_array_unsafe_get(q,2)) ))) * ( (Float)(10) )) * acc)) + 1);
HXLINE(1945)						{
HXLINE(1945)							int _g1 = 0;
HXDLIN(1945)							int _g2 = d1;
HXDLIN(1945)							while((_g1 < _g2)){
HXLINE(1945)								_g1 = (_g1 + 1);
HXDLIN(1945)								int k = (_g1 - 1);
HXLINE(1946)								Float t = ((( (Float)(k) ) / ( (Float)(d1) )) - ((Float)0.01));
HXLINE(1947)								Float s = ((( (Float)((k + 1)) ) / ( (Float)(d1) )) + ((Float)0.01));
HXLINE(1948)								Float rhs = (( (Float)(1) ) - t);
HXDLIN(1948)								Float _x = (( (Float)(_hx_array_unsafe_get(p,0)) ) * rhs);
HXDLIN(1948)								Float _y = (( (Float)(_hx_array_unsafe_get(p,1)) ) * rhs);
HXDLIN(1948)								Float _z = (( (Float)(_hx_array_unsafe_get(p,2)) ) * rhs);
HXDLIN(1948)								::Array< Float > this1 = ::Array_obj< Float >::__new(3);
HXDLIN(1948)								::Array< Float > this2 = this1;
HXDLIN(1948)								this2->__unsafe_set(0,_x);
HXDLIN(1948)								this2->__unsafe_set(1,_y);
HXDLIN(1948)								this2->__unsafe_set(2,_z);
HXDLIN(1948)								::Array< Float > this3 = this2;
HXDLIN(1948)								Float _x1 = (( (Float)(_hx_array_unsafe_get(q,0)) ) * t);
HXDLIN(1948)								Float _y1 = (( (Float)(_hx_array_unsafe_get(q,1)) ) * t);
HXDLIN(1948)								Float _z1 = (( (Float)(_hx_array_unsafe_get(q,2)) ) * t);
HXDLIN(1948)								::Array< Float > this4 = ::Array_obj< Float >::__new(3);
HXDLIN(1948)								::Array< Float > this5 = this4;
HXDLIN(1948)								this5->__unsafe_set(0,_x1);
HXDLIN(1948)								this5->__unsafe_set(1,_y1);
HXDLIN(1948)								this5->__unsafe_set(2,_z1);
HXDLIN(1948)								::Array< Float > rhs1 = this5;
HXDLIN(1948)								Float _x2 = ( (Float)(_hx_array_unsafe_get(this3,0)) );
HXDLIN(1948)								Float _x3 = (_x2 + ( (Float)(_hx_array_unsafe_get(rhs1,0)) ));
HXDLIN(1948)								Float _y2 = ( (Float)(_hx_array_unsafe_get(this3,1)) );
HXDLIN(1948)								Float _y3 = (_y2 + ( (Float)(_hx_array_unsafe_get(rhs1,1)) ));
HXDLIN(1948)								Float _z2 = ( (Float)(_hx_array_unsafe_get(this3,2)) );
HXDLIN(1948)								Float _z3 = (_z2 + ( (Float)(_hx_array_unsafe_get(rhs1,2)) ));
HXDLIN(1948)								::Array< Float > this6 = ::Array_obj< Float >::__new(3);
HXDLIN(1948)								::Array< Float > this7 = this6;
HXDLIN(1948)								this7->__unsafe_set(0,_x3);
HXDLIN(1948)								this7->__unsafe_set(1,_y3);
HXDLIN(1948)								this7->__unsafe_set(2,_z3);
HXDLIN(1948)								::Array< Float > a = this7;
HXLINE(1949)								Float rhs2 = (( (Float)(1) ) - s);
HXDLIN(1949)								Float _x4 = (( (Float)(_hx_array_unsafe_get(p,0)) ) * rhs2);
HXDLIN(1949)								Float _y4 = (( (Float)(_hx_array_unsafe_get(p,1)) ) * rhs2);
HXDLIN(1949)								Float _z4 = (( (Float)(_hx_array_unsafe_get(p,2)) ) * rhs2);
HXDLIN(1949)								::Array< Float > this8 = ::Array_obj< Float >::__new(3);
HXDLIN(1949)								::Array< Float > this9 = this8;
HXDLIN(1949)								this9->__unsafe_set(0,_x4);
HXDLIN(1949)								this9->__unsafe_set(1,_y4);
HXDLIN(1949)								this9->__unsafe_set(2,_z4);
HXDLIN(1949)								::Array< Float > this10 = this9;
HXDLIN(1949)								Float _x5 = (( (Float)(_hx_array_unsafe_get(q,0)) ) * s);
HXDLIN(1949)								Float _y5 = (( (Float)(_hx_array_unsafe_get(q,1)) ) * s);
HXDLIN(1949)								Float _z5 = (( (Float)(_hx_array_unsafe_get(q,2)) ) * s);
HXDLIN(1949)								::Array< Float > this11 = ::Array_obj< Float >::__new(3);
HXDLIN(1949)								::Array< Float > this12 = this11;
HXDLIN(1949)								this12->__unsafe_set(0,_x5);
HXDLIN(1949)								this12->__unsafe_set(1,_y5);
HXDLIN(1949)								this12->__unsafe_set(2,_z5);
HXDLIN(1949)								::Array< Float > rhs3 = this12;
HXDLIN(1949)								Float _x6 = ( (Float)(_hx_array_unsafe_get(this10,0)) );
HXDLIN(1949)								Float _x7 = (_x6 + ( (Float)(_hx_array_unsafe_get(rhs3,0)) ));
HXDLIN(1949)								Float _y6 = ( (Float)(_hx_array_unsafe_get(this10,1)) );
HXDLIN(1949)								Float _y7 = (_y6 + ( (Float)(_hx_array_unsafe_get(rhs3,1)) ));
HXDLIN(1949)								Float _z6 = ( (Float)(_hx_array_unsafe_get(this10,2)) );
HXDLIN(1949)								Float _z7 = (_z6 + ( (Float)(_hx_array_unsafe_get(rhs3,2)) ));
HXDLIN(1949)								::Array< Float > this13 = ::Array_obj< Float >::__new(3);
HXDLIN(1949)								::Array< Float > this14 = this13;
HXDLIN(1949)								this14->__unsafe_set(0,_x7);
HXDLIN(1949)								this14->__unsafe_set(1,_y7);
HXDLIN(1949)								this14->__unsafe_set(2,_z7);
HXDLIN(1949)								::Array< Float > b = this14;
HXLINE(1950)								Float o = ( (Float)(_hx_array_unsafe_get(a,2)) );
HXDLIN(1950)								Float o1 = ((o + ( (Float)(_hx_array_unsafe_get(b,2)) )) / ( (Float)(2) ));
HXLINE(1951)								int oi = ::Std_obj::_hx_int((( (Float)(255) ) - (o1 * ( (Float)(255) ))));
HXLINE(1952)								{
HXLINE(1952)									::String x = ((HX_("    <line x1=\"",2e,72,00,e3) + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(a,0)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" y1=\"",9b,09,15,55));
HXDLIN(1952)									::String x1 = ((x + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(a,1)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" x2=\"",3d,95,6c,54));
HXDLIN(1952)									::String x2 = ((x1 + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(b,0)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" y2=\"",dc,cb,15,55));
HXDLIN(1952)									::String x3 = ((((((((x2 + (( (Float)(::Math_obj::round((( (Float)(_hx_array_unsafe_get(b,1)) ) * ( (Float)(100) )))) ) / ( (Float)(100) ))) + HX_("\" stroke=\"rgb(",b6,9c,cc,2b)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(",",2c,00,00,00)) + oi) + HX_(")\"/>\n",c2,e9,0a,b2));
HXDLIN(1952)									if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1952)										out->flush();
            									}
HXDLIN(1952)									if (::hx::IsNull( out->b )) {
HXLINE(1952)										out->b = ::Array_obj< ::String >::__new(1)->init(0,::Std_obj::string(x3));
            									}
            									else {
HXLINE(1952)										::Array< ::String > out1 = out->b;
HXDLIN(1952)										out1->push(::Std_obj::string(x3));
            									}
            								}
            							}
            						}
HXLINE(1954)						gid = (gid + 1);
            					}
            				}
HXLINE(1956)				{
HXLINE(1956)					if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1956)						out->flush();
            					}
HXDLIN(1956)					if (::hx::IsNull( out->b )) {
HXLINE(1956)						out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_16,1);
            					}
            					else {
HXLINE(1956)						out->b->push(HX_("  </g>\n",a0,97,84,a9));
            					}
            				}
            			}
            		}
HXLINE(1958)		{
HXLINE(1958)			if (::hx::IsNotNull( out->charBuf )) {
HXLINE(1958)				out->flush();
            			}
HXDLIN(1958)			if (::hx::IsNull( out->b )) {
HXLINE(1958)				out->b = ::Array_obj< ::String >::fromData( _hx_array_data_04a6b36a_17,1);
            			}
            			else {
HXLINE(1958)				out->b->push(HX_("</svg>\n",1d,0c,30,0f));
            			}
            		}
HXLINE(1959)		return out->toString();
            	}


STATIC_HX_DEFINE_DYNAMIC_FUNC2(SVGWriter_obj,gradients,return )


SVGWriter_obj::SVGWriter_obj()
{
}

bool SVGWriter_obj::__GetStatic(const ::String &inName, Dynamic &outValue, ::hx::PropertyAccess inCallProp)
{
	switch(inName.length) {
	case 2:
		if (HX_FIELD_EQ(inName,"rd") ) { outValue = rd_dyn(); return true; }
		break;
	case 5:
		if (HX_FIELD_EQ(inName,"lines") ) { outValue = lines_dyn(); return true; }
		break;
	case 9:
		if (HX_FIELD_EQ(inName,"polylines") ) { outValue = polylines_dyn(); return true; }
		if (HX_FIELD_EQ(inName,"gradients") ) { outValue = gradients_dyn(); return true; }
	}
	return false;
}

#ifdef HXCPP_SCRIPTABLE
static ::hx::StorageInfo *SVGWriter_obj_sMemberStorageInfo = 0;
static ::hx::StaticInfo *SVGWriter_obj_sStaticStorageInfo = 0;
#endif

::hx::Class SVGWriter_obj::__mClass;

static ::String SVGWriter_obj_sStaticFields[] = {
	HX_("rd",b2,63,00,00),
	HX_("lines",ff,dd,01,75),
	HX_("polylines",33,0c,bc,77),
	HX_("gradients",83,78,13,cd),
	::String(null())
};

void SVGWriter_obj::__register()
{
	SVGWriter_obj _hx_dummy;
	SVGWriter_obj::_hx_vtable = *(void **)&_hx_dummy;
	::hx::Static(__mClass) = new ::hx::Class_obj();
	__mClass->mName = HX_("apparentridges.SVGWriter",6a,b3,a6,04);
	__mClass->mSuper = &super::__SGetClass();
	__mClass->mConstructEmpty = &__CreateEmpty;
	__mClass->mConstructArgs = &__Create;
	__mClass->mGetStaticField = &SVGWriter_obj::__GetStatic;
	__mClass->mSetStaticField = &::hx::Class_obj::SetNoStaticField;
	__mClass->mStatics = ::hx::Class_obj::dupFunctions(SVGWriter_obj_sStaticFields);
	__mClass->mMembers = ::hx::Class_obj::dupFunctions(0 /* sMemberFields */);
	__mClass->mCanCast = ::hx::TCanCast< SVGWriter_obj >;
#ifdef HXCPP_SCRIPTABLE
	__mClass->mMemberStorageInfo = SVGWriter_obj_sMemberStorageInfo;
#endif
#ifdef HXCPP_SCRIPTABLE
	__mClass->mStaticStorageInfo = SVGWriter_obj_sStaticStorageInfo;
#endif
	::hx::_hx_RegisterClass(__mClass->mName, __mClass);
}

} // end namespace apparentridges
