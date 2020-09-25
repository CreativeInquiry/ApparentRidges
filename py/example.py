from apparentridges import *

mesh = apparentridges_OBJParser.fromFile("../testdata/lucy.obj");
render = apparentridges_Render(mesh,800,800);
render.autoPlace();
render.apparentRidges(0.1);
render.buildPolylines();

rstr = apparentridges_SVGWriter.polylines(render);
open("out.svg",'w').write(rstr);