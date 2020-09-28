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
**/obj/**/*.*
**/obj/*.*
hl/**/hl/*.*
**/_*/**/*.*
**/_*.*\
""".split("\n");

def realglob(ignores):
  G = []
  for ig in ignores:
    for i in range(0,5):
      iig = ig.replace("**/","*/"*i,1);
      for j in range(0,5):
        ijg = iig.replace("**/","*/"*j);
        ggg = glob(ijg);
        G += ggg
        # print(ig,ijg,ggg);
  return G

nof = list(set(realglob(ignores)))

open(".gitattributes",'w').write(
  "\n".join([x+" linguist-generated\n"+x+" linguist-detectable=false" for x in nof])
)

allf = """\
*.*
**/*.*
**/**/*.*
**/**/**/*.*
**/**/**/**/*.*\
""".split("\n");

gigf = open(".gitignore",'r').read().split("\n");
gigf = gigf + [x+"/**/*.*" for x in gigf]
gigf = list(set(realglob(gigf)))

yesf = list(set(realglob(allf)))
yesf = [x for x in yesf if x not in nof]
yesf = [x for x in yesf if x not in gigf]
yesf = [x for x in yesf if "docs" not in x]
yesf = [x for x in yesf if "svg" not in x]
yesf = [x for x in yesf if "xml" not in x]

open(".gitattributes",'w').write(
   "\n".join([x+" linguist-generated\n"+x+" linguist-detectable=false" for x in nof])
  +"\n".join([x+" linguist-detectable=true" for x in yesf])
  +"\n*.java linguist-language=Java"
  +"\n*.c linguist-detectable=false"
  +"\n*.h linguist-detectable=false"
)