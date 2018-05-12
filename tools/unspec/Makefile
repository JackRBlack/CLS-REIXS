help:
	@echo This Makefile compiles program unspec
	@echo Usage:
	@echo	make g++     ... g++ under unix or unixish environment
	@echo	make static  ... g++ as above, static binary
	@echo	make clang++ ... clang++ compiler
	@echo	make emx     ... gcc/emx under OS/2 or DOS
	@echo	make hpcc    ... HP CC under HP-UX
	@echo	make deccxx  ... DEC CXX under VMS

### Makefile body:

g++:
	g++ -Wno-deprecated -O2 -s unspec.cpp -o unspec

static:
	g++ -O2 -s unspec.cpp -o unspec -static

clang++:
	clang++ -O2 -s unspec.cpp -o unspec

emx:
	gcc -Wall -O2 -s unspec.cpp -o unspec.exe -lstdcpp

hpcc:
	CC -O -s -w unspec.cpp -o unspec

deccxx:
	ren unspec.cpp unspec.cxx
	CXX unspec.cxx
	LINK /exe=unspec.exe unspec.obj

debug:
	g++ -Wno-deprecated -Wall -g unspec.cpp -o unspec

