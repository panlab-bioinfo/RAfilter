object = main.o kmerfind4.o filter.o
target = rafilter
LIBS += -L htslib/ -lhts -I htslib/htslib/


rafilter : $(object)
	echo Generate an executable file
	g++ -Wl,rpath=htslib/ -o $(target) $(object) $(LIBS) -lpthread

filter.o:
	cd htslib && autoreconf -i
	cd htslib && ./configure
	cd htslib && make
	g++ -c filter.cpp $(LIBS)

.PHONY : clean
clean :
	cd htslib && make clean
	rm $(object) $(target)
