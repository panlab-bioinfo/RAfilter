object = main.o kmerfind4.o filter.o
target = rafilter
CXXFLAGS=-Wl,--rpath=$(PWD)/htslib/ --std=c++11 -L htslib/ -lhts -static-libstdc++ -I htslib/htslib/ -lpthread

$(target) : $(object)
	@echo Generate an executable file
	g++ -o $(target) $(object) $(CXXFLAGS)

filter.o:
	cd htslib && autoreconf -i
	cd htslib && ./configure --disable-bz2
	cd htslib && make
	g++ -c filter.cpp $(CXXFLAGS)

info :
	@echo The program "rafilter" will obtained in the processing.
	@echo pwd is $(PWD)
	@echo compiler is $(CC)
	@echo  Flags is $(CXXFLAGS)

.PHONY : clean
clean :
	cd htslib && make clean
	rm $(object) $(target)
