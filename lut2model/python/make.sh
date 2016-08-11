swig -c++ -python lut2model.i 
NUMPYFLAGS=-I/home/savchenk/.pythonbrew/pythons/Python-2.7.3/lib/python2.7/site-packages/numpy/core/include/
#NUMPYFLAGS=-I/usr/local/python/python-2.7/lib/python2.7/site-packages/numpy/core/include
g++ -fPIC -pthread -fno-strict-aliasing -O3 -DNDEBUG -fwrapv -Wall -Wstrict-prototypes `python-config --cflags` $NUMPYFLAGS -c lut2model_wrap.cxx -o lut2model_wrap.o
g++ -fPIC -pthread -fno-strict-aliasing -O3 -DNDEBUG -fwrapv -Wall `python-config --cflags` $NUMPYFLAGS -c lut2model.cxx -o lut2model.o
g++ -pthread -O3 -shared lut2model_wrap.o lut2model.o `python-config --ldflags` -L/home/savchenk/.pythonbrew/pythons/Python-2.7.3/lib -o _lut2model.so
