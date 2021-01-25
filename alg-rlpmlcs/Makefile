OBJS += rlpmlcs.o

librlpmlcs.a: $(OBJS)
	${AR} -r -o $@ $^
	cp librlpmlcs.a ${LIBSDIR}

%.o: %.cpp
	${CC} $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a