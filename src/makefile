object = main.o kmerfind4.o filter.o
target = rafilter
CXXFLAGS=-Wl,--rpath=$(PWD)/htslib/  -L htslib/ -lhts -I htslib/htslib/ -lpthread

$(target) : $(object)
	@echo Generate an executable file
	g++ -o $(target) $(object) $(CXXFLAGS)

filter.o:
	cd htslib && autoreconf -i
	cd htslib && ./configure
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
