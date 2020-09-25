all: libs clis tests

libs: java js lua py cpp

clis: cli-py cli-cpp cli-hlc

_:

java: _
	haxe apparentridges/ApparentRidges.hx -D real_position --java java/ApparentRidges

js: _
	haxe apparentridges/ApparentRidges.hx --js js/apparentridges.js

lua: _
	haxe apparentridges/ApparentRidges.hx -D lua_jit -D lua_vanilla --lua lua/apparentridges.lua

py: _
	haxe apparentridges/ApparentRidges.hx --python py/apparentridges.py

c: _
	haxe apparentridges/ApparentRidges.hx -hl hl/c/main.c
	cd hl/c; gcc -O3 main.c -o out.o -I. -lhl -I/usr/local/Cellar/hashlink/1.11_1/libexec/include -L /usr/local/Cellar/hashlink/1.11_1/libexec/lib;
	cd hl/c; ar rcs libapparentridges.a out.o

cpp: _
	haxe apparentridges/ApparentRidges.hx --cpp cpp/ApparentRidges/ -D static_link

cpp-example: _
	cd cpp; g++ example.cpp -I/usr/local/lib/haxe/lib/hxcpp/4,1,15/include -I ./ApparentRidges/include --std=c++11 -L./ApparentRidges/ -loutput

cli-py: _
	haxe --main Cli.hx --python py/cli.py

cli-cpp: _
	haxe --main Cli.hx --cpp cpp/cli/

cli-hlc: _
	haxe -hl hl/cli/main.c -main Cli;
	cd hl/cli; gcc -O3 main.c -o cli -I. -lhl -I/usr/local/Cellar/hashlink/1.11_1/libexec/include -L /usr/local/Cellar/hashlink/1.11_1/libexec/lib;

tests: _
	node js/generatesamples.js all

test-%: _
	node js/generatesamples.js $*