local ar = require "apparentridges"
local io = require "io"

local mesh = ar.apparentridges.OBJParser.fromFile("../testdata/lucy.obj");
local render = ar.apparentridges.Render.new(mesh,800,800);
render:autoPlace();
render:apparentRidges(0.1);
render:buildPolylines();

local rstr = ar.apparentridges.SVGWriter.polylines(render);
local file = io.open("out.svg", "w")
io.output(file)
io.write(rstr)
io.close(file)