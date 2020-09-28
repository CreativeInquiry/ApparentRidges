from glob import glob
ignores="""\
**/haxe/**/*.*
**/haxe/*.*
**/types/**/*.*
**/types/*.*
**/std/**/*.*
**/std/*.*
**/_std/**/*.*
**/_std/*.*
**/sys/**/*.*
**/sys/*.*
hl/**/hl/*.*
**/_*/**/*.*
**/_*.*
""".split("\n");

G = []
for ig in ignores:
  for i in range(0,5):
    iig = ig.replace("**/","*/"*i,1);
    for j in range(0,5):
      ijg = iig.replace("**/","*/"*j);
      ggg = glob(ijg);
      G += ggg
      # print(ig,ijg,ggg);
open(".gitattributes",'w').write(
  "\n".join([x+" linguist-generated\n"+x+" linguist-detectable=false" for x in G])
)