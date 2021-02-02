CPPFLAGS += -std=c++11
CPPFLAGS += -fpermissive
CPPFLAGS += -I ${shell pwd}/include/

CC = g++
AR = ar
LD = ld

.DEFAULT_GOAL := all

OBJS += common/io.o common/phash.o 
OBJS += main.o
OBJS += utils/tool.o

LIBS := 
LIBSFLAGS :=
LIBSDIR := ${shell pwd}/libs/

# alg program tag start
LIBS += alg-topmlcs/libtopmlcs.a
LIBSFLAGS += -ltopmlcs
alg-topmlcs/libtopmlcs.a:
	cd alg-topmlcs; make; cd ..

topmlcs:
	cd alg-topmlcs; make; cd ..

LIBS += alg-quickdp/libquickdp.a
LIBSFLAGS += -lquickdp
alg-quickdp/libquickdp.a:
	cd alg-quickdp; make; cd ..

quickdp:
	cd alg-quickdp; make; cd ..

LIBS += alg-promlcs/libpromlcs.a
LIBSFLAGS += -lpromlcs
alg-promlcs/libpromlcs.a:
	cd alg-promlcs; make; cd ..

promlcs:
	cd alg-promlcs; make; cd ..

LIBS += alg-mlcsapp/libmlcsapp.a
LIBSFLAGS += -lmlcsapp
alg-mlcsapp/libmlcsapp.a:
	cd alg-mlcsapp; make; cd ..

mlcsapp:
	cd alg-mlcsapp; make; cd ..

LIBS += alg-hasmlcs/libhasmlcs.a
LIBSFLAGS += -lhasmlcs
alg-hasmlcs/libhasmlcs.a:
	cd alg-hasmlcs; make; cd ..

hasmlcs:
	cd alg-hasmlcs; make; cd ..

# alg program tag end

TARGET = LCSCalculator

export CC CPPFLAGS LD AR LIBS LIBSDIR

$(TARGET): $(OBJS) $(LIBS)
	$(CC) -o $(TARGET) $(OBJS) -L$(LIBSDIR) ${LIBSFLAGS}

%.o : %.cpp
	$(CC) $(CPPFLAGS) -o $@ $< -c

genTest: utils/genTest.o
	$(CC) -o genTest utils/genTest.o $(CPPFLAGS)

all: genTest $(TARGET)

update-lib:
# end target update-lib

clean:
	for i in `ls` ; do \
		if [[ $$i == alg-* ]]; then \
    		cd $$i && make clean && cd .. ; \
		fi ; \
	done
	cd common; make clean; cd ..
	cd utils; make clean; cd ..
	cd libs; rm *.a; cd ..
	rm *.o $(TARGET) genTest
