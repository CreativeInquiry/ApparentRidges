import apparentridges.ApparentRidges;

class Cli{
  
  static public function main():Int {
    inline function shortId() : String{
      var id : String = "";
      for (i in 0...4){
        id += String.fromCharCode(Std.int(Math.random()*26)+0x61);
      }
      return id;
    }

                                                                                                                    
    
    var abbrev : Map<String, String> = ["-o"=>"--output","-w"=>"--width","-h"=>"--height"];
    var argv : Array<String> = Sys.args();
    var kwargs:Map<String, String> = [
      "--width"=>"800",
      "--height"=>"800",
      "--output"=>("out-"+shortId()+".svg"),
      "--transform"=>"auto()",
      "--cull"=>"true",
      "--thresh"=>"0.1",
      "--mode"=>"gradients",
      "--verbose"=>"1",
    ];
    var desc:Map<String, String> = [
      "--width"=>"width of output image",
      "--height"=>"height of output image",
      "--output"=>"output filename\n(randomly generated if not specified)",
      "--transform"=>"place the model before rendering\n(multiple commands are left-multiplied in order)\navailable commands:\nfocal(100) translate(1,2,3) scale(4,4,4)\nrotateX(1) rotateY(2) rotateZ(3);\nmatrix(11,12,13,14,21,22,...,43,44)\nrotation angles are in degrees\nuse auto() or auto(zFactor,fFactor) for automatic\nplacement. (cam is fixed at (0,0,0) pointing at +Z,\nfocal() needs to be specified for manual placement)",
      "--cull"=>"don't draw occluded faces\noptions: true/false/custom float value (e.g. 1.5)\nthe float being a multiplier to the 'epsilon',\nallowing ridges just barely occluded to show up",
      "--thresh"=>"apparent ridge threshold\nsmaller => more detailed, larger => cleaner",
      "--mode"=>"visualization technique, options:\nvertices, edges, \nlines (disconnected ridge segments),\npolylines (connected ridges via post-proc),\ngradients (nice render showing ridge strength),\npixels (raster out: depth, norm., curv. maps)",
      "--verbose"=>"verbosity: 0 => errors only, 1 => logs",
    ];

    if (argv.length == 0 || argv[0] == "--help"){
      Sys.println("
       ___           ___      
      /\\  \\         /\\  \\     
     /  \\  \\       /  \\  \\    
    / /\\ \\  \\     / /\\ \\  \\   
   /  \\~\\ \\  \\   /  \\~\\ \\  \\  
  / /\\ \\ \\ \\__\\ / /\\ \\ \\ \\__\\ 
  \\/__\\ \\/ /  / \\/_|  \\/ /  / 
       \\  /  /     | |  /  /  
       / /  /      | |\\/__/   
      / /  /       | |  |     
      \\/__/ PPARENT \\|__| IDGES (CLI)

Render line drawings of 3D meshes.

usage:   apparentridges [options] [.obj files...]

option        default       description
");
      for (k in kwargs.keys()){
        var defau : String = kwargs[k];
        var ln : String = k;
        for (a in abbrev.keys()){
          if (abbrev[a] == k){
            ln+=", "+a;
          }
        }
        ln = StringTools.rpad(ln," ",14);
        ln += defau;
        ln = StringTools.rpad(ln," ",28);
        if (desc.exists(k)){
          ln += StringTools.replace(desc[k],"\n","\n"+([for (i in 0...28) " "].join("")));
        }
        ln += "\n";
        Sys.println(ln);
      }
      
      return 0;
    }

    var args:Array<String> = [];

    var i : Int = 0;
    while (i < argv.length){
      if (argv[i].length == 0){
        i++;
        continue;
      }
      if (argv[i].charAt(0)=="-"){
        if (i+1 >= argv.length){
          trace("[error] option without value: "+argv[i]);
          return 1;
        }
        var key : String = argv[i];
        var val : String = argv[i+1];
        if (abbrev.exists(key)){
          key = abbrev[key];
        }
        kwargs[key] = StringTools.replace(val," ","");
        i+=2;
        continue;
      }
      args.push(argv[i]);
      i++;
    }
    if (args.length < 1){
      trace("[warn] no input file.");
      return 0;
    }

    inline function floatFuncArgs(s : String) : Array<Float> {
      return [for (x in s.split("(")[1].split(")")[0].split(",")) Std.parseFloat(x)].filter(function(x : Float):Bool{
        return !Math.isNaN(x);
      });
      
    }

    var verb : Int = Std.parseInt(kwargs["--verbose"]);
    
    for (inp in args){
      if (verb>0)trace("[info] working on "+inp+"...");
      var mesh : Mesh = OBJParser.fromFile(inp);
      var render : Render = new Render(
        mesh,Std.parseInt(kwargs["--width"]),Std.parseInt(kwargs["--height"]));

      render.setVerbose(verb);

      var m : Array<Float> = Util.matIden();
      for (t in kwargs["--transform"].split(")")){
        if (t.length==0){
          continue;
        }
        var nums : Array<Float> = floatFuncArgs(t);

        if (StringTools.startsWith(t,"rotateX")){
          m = Util.matMult(Util.matRotx(nums[0]*Math.PI/180), m);
        }else if (StringTools.startsWith(t,"rotateY")){
          m = Util.matMult(Util.matRoty(nums[0]*Math.PI/180), m);
        }else if (StringTools.startsWith(t,"rotateZ")){
          m = Util.matMult(Util.matRotz(nums[0]*Math.PI/180), m);
        }else if (StringTools.startsWith(t,"scale")){
          m = Util.matMult(Util.matScal(nums[0],nums[1],nums[2]), m);
        }else if (StringTools.startsWith(t,"translate")){
          m = Util.matMult(Util.matTrsl(nums[0],nums[1],nums[2]), m);
        }else if (StringTools.startsWith(t,"matrix")){
          m = Util.matMult(nums, m);
        }else if (StringTools.startsWith(t,"focal")){
          render.setFocal(nums[0]);
        }else if (StringTools.startsWith(t, "auto")){
          render.transform(m);
          m = Util.matIden();
          if (nums.length == 0){
            render.autoPlace();
          }else{
            render.autoPlace(nums[0],nums[1]);
          }
        }else{
          trace("[warn] invalid transform expression: "+t);
        }
      }
      render.transform(m);
      
      
      if (kwargs["--mode"] != "edges" && kwargs["--mode"] != "vertices"){
        if (kwargs["--cull"]=="true"){
          render.apparentRidges(Std.parseFloat(kwargs["--thresh"]));
        }else if (kwargs["--cull"]=="false"){
          render.apparentRidges(Std.parseFloat(kwargs["--thresh"]),-1);
        }else{
          render.apparentRidges(Std.parseFloat(kwargs["--thresh"]),Std.parseFloat(kwargs["--cull"]));
        }
      }

      if (verb>0)trace("[info] generating output...");

      var out : String = "";

      if (kwargs["--mode"]=="edges"){
        render.edges();
        out = SVGWriter.lines(render);
      }else if (kwargs["--mode"]=="vertices"){
        render.vertices();
        out = SVGWriter.lines(render);
      }else if (kwargs["--mode"]=="lines"){
        out = SVGWriter.lines(render);
      }else if (kwargs["--mode"]=="polylines"){
        render.buildPolylines();
        out = SVGWriter.polylines(render);
      }else if (kwargs["--mode"]=="gradients"){
        render.buildPolylines();
        out = SVGWriter.gradients(render);
      }else if (kwargs['--mode']=="pixels"){
        var d : haxe.ds.Vector<Float>;
        var o : String;
        var p : String;
        if (verb>0)trace("[info] making depth map...");
        d = PixelMap.depth(render,true);
        if (verb>0)trace("[info] making ppm string...");
        o = PixelMap.toPPMString(d,render.width,render.height,0,1);
        p = inp+".depth.ppm";
        sys.io.File.saveContent(p,o);
        if (verb>0)trace("[info] wrote "+p);

        if (verb>0)trace("[info] making normal map...");
        d = PixelMap.normal(render);
        o = PixelMap.toPPMString(d,render.width,render.height,-1,1);
        p = inp+".normal.ppm";
        sys.io.File.saveContent(p,o);
        if (verb>0)trace("[info] wrote "+p);

        if (verb>0)trace("[info] making curvature map...");
        d = PixelMap.curvature(render);
        o = PixelMap.toPPMString(d,render.width,render.height,-0.5,0.5);
        p = inp+".curv.ppm";
        sys.io.File.saveContent(p,o);
        if (verb>0)trace("[info] wrote "+p);
        
        continue;

      }else{
        trace("[error] invalid mode");
        return 1;
      }
      if (verb>0)trace("[info] writing file...");

      var path : String = kwargs["--output"];
      if (path == "stdout"){
        Sys.println(out);
      }else{
        if (args.length>1){
          path = inp+".svg";
        }
        sys.io.File.saveContent(path,out);
        if (verb>0)trace("[info] wrote "+path);
      }
    }

    if (verb>0)trace("[info] all files processed.");
    return 0;
  }
}