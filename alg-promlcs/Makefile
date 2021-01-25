OBJS += dtree.o promlcs.o

libpromlcs.a: $(OBJS)
	${AR} -r -o $@ $^
	cp libpromlcs.a ${LIBSDIR}

%.o: %.cpp
	${CC} $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a