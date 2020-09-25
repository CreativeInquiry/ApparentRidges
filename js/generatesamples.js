const fs = require('fs');
const {apparentridges} = require("./apparentridges");
const AR = apparentridges;

const VIS_MODE = "gradients"

global.planck = function(){
  let objstr = fs.readFileSync("./testdata/planck_100K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRotx(-Math.PI*0.5));
  render.transform(AR.Util.matRoty(Math.PI*0.2));
  render.transform(AR.Util.matRotx(-0.3));
  render.autoPlace(3,2.4);
  render.apparentRidges(0.25,10);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/planck.svg",rstr);
}

global.dragon = function(){
  let objstr = fs.readFileSync("./testdata/dragon_100K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRoty(Math.PI*1.0));
  render.transform(AR.Util.matRotx(-0.4));

  render.autoPlace();
  render.apparentRidges(0.08);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/dragon.svg",rstr);
}

global.buddha = function(){
  let objstr = fs.readFileSync("./testdata/buddha_100K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRoty(Math.PI*1.0));

  render.autoPlace();
  render.apparentRidges(0.1);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/buddha.svg",rstr);
}


global.bunny = function(){
  let objstr = fs.readFileSync("./testdata/bunny_70K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRoty(Math.PI*1.05));
  render.transform(AR.Util.matRotx(-0.25));
  render.autoPlace(2,1.5);

  render.apparentRidges(0.15);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/bunny.svg",rstr);
}


global.lucy = function(){
  let objstr = fs.readFileSync("./testdata/lucy_500K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRotx(-Math.PI/2));
  render.autoPlace();

  render.apparentRidges(0.4);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/lucy.svg",rstr);
}

global.nefertiti = function(){
  let objstr = fs.readFileSync("./testdata/nefertiti_100K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRotx(-Math.PI*0.5));
  render.transform(AR.Util.matRoty(Math.PI*1.1));
  render.autoPlace(3,2.5);

  render.apparentRidges(0.1);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/nefertiti.svg",rstr);

}

global.ogre = function(){
  let objstr = fs.readFileSync("./testdata/ogre_40K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);
  render.transform(AR.Util.matRoty(Math.PI*1.1));
  render.autoPlace(3,2.5);

  render.apparentRidges(0.05,5);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/ogre.svg",rstr);

}
global.armadillo = function(){
  let objstr = fs.readFileSync("./testdata/armadillo_300K.obj").toString();
  let mesh = AR.OBJParser.fromString(objstr);
  let render = new AR.Render(mesh,800,800);

  render.autoPlace();

  render.apparentRidges(0.05);
  render.buildPolylines();
  let rstr = AR.SVGWriter[VIS_MODE](render);
  fs.writeFileSync("output/armadillo.svg",rstr);

}

global.all = function(){
  console.log("drawing lucy...");
  lucy();
  console.log("drawing bunny...");
  bunny();
  console.log("drawing dragon...");
  dragon();
  console.log("drawing buddha...");
  buddha();
  console.log("drawing planck...");
  planck();
  console.log("drawing nefertiti...");
  nefertiti();
  console.log("drawing ogre...");
  ogre();
  console.log("drawing armadillo...");
  armadillo();
}

global[process.argv[2]]();