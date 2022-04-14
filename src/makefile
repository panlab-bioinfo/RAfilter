object = main.o kmerfind4.o filter.o
target = kmfilter
kmfilter: $(object)
	g++ -o $(target) $(object) -lpthread
	mv $(target) ../
.PHONY : clean
clean :
	rm $(objects)
