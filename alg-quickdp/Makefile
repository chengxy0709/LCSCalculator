OBJS += quickdp.o

libquickdp.a: $(OBJS)
	${AR} -r -o $@ $^
	cp libquickdp.a ${LIBSDIR}

%.o: %.cpp
	${CC} $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a