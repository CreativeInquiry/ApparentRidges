// Generated by Haxe 4.1.3
#ifndef INCLUDED_apparentridges__ApparentRidges_Face_Impl_
#define INCLUDED_apparentridges__ApparentRidges_Face_Impl_

#ifndef HXCPP_H
#include <hxcpp.h>
#endif

HX_DECLARE_CLASS2(apparentridges,_ApparentRidges,Face_Impl_)

namespace apparentridges{
namespace _ApparentRidges{


class HXCPP_CLASS_ATTRIBUTES Face_Impl__obj : public ::hx::Object
{
	public:
		typedef ::hx::Object super;
		typedef Face_Impl__obj OBJ_;
		Face_Impl__obj();

	public:
		enum { _hx_ClassId = 0x2f5499b8 };

		void __construct();
		inline void *operator new(size_t inSize, bool inContainer=false,const char *inName="apparentridges._ApparentRidges.Face_Impl_")
			{ return ::hx::Object::operator new(inSize,inContainer,inName); }
		inline void *operator new(size_t inSize, int extra)
			{ return ::hx::Object::operator new(inSize+extra,false,"apparentridges._ApparentRidges.Face_Impl_"); }

		inline static ::hx::ObjectPtr< Face_Impl__obj > __new() {
			::hx::ObjectPtr< Face_Impl__obj > __this = new Face_Impl__obj();
			__this->__construct();
			return __this;
		}

		inline static ::hx::ObjectPtr< Face_Impl__obj > __alloc(::hx::Ctx *_hx_ctx) {
			Face_Impl__obj *__this = (Face_Impl__obj*)(::hx::Ctx::alloc(_hx_ctx, sizeof(Face_Impl__obj), false, "apparentridges._ApparentRidges.Face_Impl_"));
			*(void **)__this = Face_Impl__obj::_hx_vtable;
			return __this;
		}

		static void * _hx_vtable;
		static Dynamic __CreateEmpty();
		static Dynamic __Create(::hx::DynamicArray inArgs);
		//~Face_Impl__obj();

		HX_DO_RTTI_ALL;
		static bool __GetStatic(const ::String &inString, Dynamic &outValue, ::hx::PropertyAccess inCallProp);
		static void __register();
		bool _hx_isInstanceOf(int inClassId);
		::String __ToString() const { return HX_("Face_Impl_",1d,c9,d8,54); }

		static ::Array< int > _new(int a,int b,int c);
		static ::Dynamic _new_dyn();

		static int get(::Array< int > this1,int i);
		static ::Dynamic get_dyn();

		static int set(::Array< int > this1,int i,int v);
		static ::Dynamic set_dyn();

		static int indexOf(::Array< int > this1,int v);
		static ::Dynamic indexOf_dyn();

};

} // end namespace apparentridges
} // end namespace _ApparentRidges

#endif /* INCLUDED_apparentridges__ApparentRidges_Face_Impl_ */ 
