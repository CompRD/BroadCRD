PYTHON_INCLUDE=/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.6.4/include/python2.6
SWIG=~ghall/bin/swig

rm -f scaffold.py _scaffold.so scaffold_wrap.o scaffold_wrap.cxx

$SWIG -c++ -I.. -I/broad/software/free/Linux/redhat_5_x86_64/pkgs/gcc_4.4.3/include/c++/4.4.3 -python scaffold.i

g++ -fPIC -Wno-deprecated -c scaffold_wrap.cxx -I.. -I$PYTHON_INCLUDE

g++ scaffold_wrap.o scaffold.a -shared -o _scaffold.so
