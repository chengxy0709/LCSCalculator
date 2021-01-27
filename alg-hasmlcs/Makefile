OBJS += hasmlcs.o priority.o

libhasmlcs.a: $(OBJS)
	$(AR) -r -o $@ $^
	cp libhasmlcs.a ${LIBSDIR}

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a