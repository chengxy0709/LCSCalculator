CPPFLAGS += -std=c++11
CPPFLAGS += -fpermissive
CPPFLAGS += -I ./include/

CC = g++
OBJS += common/io.o common/phash.o 
OBJS += mlcsapp/MLCSAPP.o
OBJS += promlcs/PRO_MLCS.o promlcs/dtree.o
OBJS += quickdp/QuickDP.o
OBJS += rlpmlcs/RLP_MLCS.o
OBJS += topmlcs/TOP_MLCS.o
OBJS += wdag/WDAG.o
OBJS += utils/tool.o
OBJS += main.o
TARGET = LCSCalculator

export CC CPPFLAGS

%.o : %.cpp
	$(CC) $(CPPFLAGS) -o $@ $< -c

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CPPFLAGS)

genTest: utils/genTest.o
	$(CC) -o genTest genTest.o $(CPPFLAGS)

clean:
	cd common; rm *.o; cd ..
	cd mlcsapp; rm *.o; cd ..
	cd promlcs; rm *.o; cd ..
	cd quickdp; rm *.o; cd ..
	cd rlpmlcs; rm *.o; cd ..
	cd topmlcs; rm *.o; cd ..
	cd utils; rm *.o; cd ..
	cd wdag; rm *.o; cd ..
	rm *.o $(TARGET) genTest