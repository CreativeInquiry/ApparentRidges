// Generated by Haxe 4.1.3
#ifndef INCLUDED_apparentridges_BBox
#define INCLUDED_apparentridges_BBox

#ifndef HXCPP_H
#include <hxcpp.h>
#endif

HX_DECLARE_CLASS1(apparentridges,BBox)

namespace apparentridges{


class HXCPP_CLASS_ATTRIBUTES BBox_obj : public ::hx::Object
{
	public:
		typedef ::hx::Object super;
		typedef BBox_obj OBJ_;
		BBox_obj();

	public:
		enum { _hx_ClassId = 0x4272c92a };

		void __construct();
		inline void *operator new(size_t inSize, bool inContainer=true,const char *inName="apparentridges.BBox")
			{ return ::hx::Object::operator new(inSize,inContainer,inName); }
		inline void *operator new(size_t inSize, int extra)
			{ return ::hx::Object::operator new(inSize+extra,true,"apparentridges.BBox"); }
		static ::hx::ObjectPtr< BBox_obj > __new();
		static ::hx::ObjectPtr< BBox_obj > __alloc(::hx::Ctx *_hx_ctx);
		static void * _hx_vtable;
		static Dynamic __CreateEmpty();
		static Dynamic __Create(::hx::DynamicArray inArgs);
		//~BBox_obj();

		HX_DO_RTTI_ALL;
		::hx::Val __Field(const ::String &inString, ::hx::PropertyAccess inCallProp);
		::hx::Val __SetField(const ::String &inString,const ::hx::Val &inValue, ::hx::PropertyAccess inCallProp);
		void __GetFields(Array< ::String> &outFields);
		static void __register();
		void __Mark(HX_MARK_PARAMS);
		void __Visit(HX_VISIT_PARAMS);
		bool _hx_isInstanceOf(int inClassId);
		::String __ToString() const { return HX_("BBox",e9,8a,d2,2b); }

		::Array< Float > min;
		::Array< Float > max;
		::Array< Float > centroid();
		::Dynamic centroid_dyn();

		void add(::Array< Float > p);
		::Dynamic add_dyn();

		void merge( ::apparentridges::BBox bb);
		::Dynamic merge_dyn();

		Float surfaceArea();
		::Dynamic surfaceArea_dyn();

};

} // end namespace apparentridges

#endif /* INCLUDED_apparentridges_BBox */ 
