// g++ example.cpp -I/usr/local/lib/haxe/lib/hxcpp/4,1,15/include -I ./ApparentRidges/include --std=c++11 -L./ApparentRidges/ -loutput

#include "ApparentRidges/include/apparentridges/Mesh.h"
#include "ApparentRidges/include/apparentridges/Render.h"
#include "ApparentRidges/include/apparentridges/OBJParser.h"
#include "ApparentRidges/include/apparentridges/SVGWriter.h"
#include "ApparentRidges/include/apparentridges/Util.h"

int main(){
  HX_TOP_OF_STACK
  ::hx::Boot();
	__boot_all();
  hx::StackContext *_hx_ctx = hx::gMainThreadContext;

  ::apparentridges::Mesh mesh = ::apparentridges::OBJParser_obj::fromFile(HX_("../testdata/lucy_100K.obj",ff,ff,ff,ff));
  ::apparentridges::Render render =  ::apparentridges::Render_obj::__alloc(HX_CTX,mesh,800,800);
  render->autoPlace(null(),null());
  render->apparentRidges(0.1,null());
  String s = ::apparentridges::SVGWriter_obj::lines(render,true);
  ::apparentridges::Util_obj::writeFile(HX_("out.svg",ff,ff,ff,ff),s);

  return 0;
}
