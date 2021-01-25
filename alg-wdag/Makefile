OBJS += wdag.o

libwdag.a: $(OBJS)
	${AR} -r -o $@ $^
	cp libwdag.a ${LIBSDIR}

%.o: %.cpp
	${CC} $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o