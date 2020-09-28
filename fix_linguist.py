from glob import glob
ignores="""\
**/haxe/**/*.*
**/haxe/*.*
**/std/**/*.*
**/std/*.*
**/_std/**/*.*
**/_std/*.*
**/sys/**/*.*
**/sys/*.*
hl/**/hl/*.*
""".split("\n");

G = []
for ig in ignores:
  for i in range(0,4):
    iig = ig.replace("**/","*/"*i,1);
    for j in range(0,4):
      ijg = iig.replace("**/","*/"*j);
      ggg = glob(ijg);
      G += ggg
      # print(ig,ijg,ggg);
open(".gitattributes",'w').write(
  "\n".join([x+" linguist-detectable=false" for x in G])
)