object = main.o kmerfind4.o filter.o
target = kmfilter
LIBS += -L htslib/ -lhts -I htslib/htslib/


kmfilter : $(object)
	echo Generate an executable file
	g++ -g -o $(target) $(object) $(LIBS) -lpthread
	echo add htslib to LD_LIBRARY_PATH
	echo export LD_LIBRARY_PATH=`pwd`/htslib/ >> ~/.bashrc

filter.o:
	cd htslib && autoreconf -i
	cd htslib && ./configure
	cd htslib && make
	g++ -c filter.cpp $(LIBS)

.PHONY : clean
clean :
	cd htslib && make clean
	rm $(object) $(target)
