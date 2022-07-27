object = main.o kmerfind4.o filter.o
target = kmfilter
LIBS += -L  htslib/ -lhts
kmfilter: $(object)
	g++ -g -o $(target) $(object) $(LIBS) -lpthread


.PHONY : clean
clean :
	rm $(objects)
