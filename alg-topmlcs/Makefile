OBJS += topmlcs.o

libtopmlcs.a: $(OBJS)
	${AR} -r -o $@ $^
	cp libtopmlcs.a ${LIBSDIR}

%.o: %.cpp
	${CC} $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a