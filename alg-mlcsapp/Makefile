OBJS += mlcsapp.o

libmlcsapp.a: $(OBJS)
	$(AR) -r -o $@ $^
	cp libmlcsapp.a ${LIBSDIR}

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a