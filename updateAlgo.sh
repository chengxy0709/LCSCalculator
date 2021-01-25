#!/bin/bash

ALG_OBJ_S=`grep "alg program tag start" -n ./Makefile  | awk -F: '{print $1}'`
ALG_OBJ_E=`grep "alg program tag end" -n ./Makefile  | awk -F: '{print $1}'`
DEL_ALG_OBJ_S=`expr $ALG_OBJ_S + 1`
DEL_ALG_OBJ_E=`expr $ALG_OBJ_E - 1`

ALG_HEADER_S=`grep "alg program header start" -n ./include/io.h  | awk -F: '{print $1}'`
ALG_HEADER_E=`grep "alg program header end" -n ./include/io.h  | awk -F: '{print $1}'`
DEL_ALG_HEADER_S=`expr $ALG_HEADER_S + 1`
DEL_ALG_HEADER_E=`expr $ALG_HEADER_E - 1`

ALG_INTERFACE_S=`grep "alg program interface start" -n ./common/io.cpp  | awk -F: '{print $1}'`
ALG_INTERFACE_E=`grep "alg program interface end" -n ./common/io.cpp  | awk -F: '{print $1}'`
DEL_ALG_INTERFACE_S=`expr $ALG_INTERFACE_S + 1`
DEL_ALG_INTERFACE_E=`expr $ALG_INTERFACE_E - 1`

# delete all alg program objects in Makefile
if [ $DEL_ALG_OBJ_E -gt $DEL_ALG_OBJ_S ]; then
    sed -i "$DEL_ALG_OBJ_S,$DEL_ALG_OBJ_E d" Makefile 
fi

# delete all alg program headers in io.h
if [ $DEL_ALG_HEADER_E -gt $DEL_ALG_HEADER_S ]; then
    sed -i "$DEL_ALG_HEADER_S,$DEL_ALG_HEADER_E d" ./include/io.h
fi

# delete all alg program interfaces in io.cpp
if [ $DEL_ALG_INTERFACE_E -gt $DEL_ALG_INTERFACE_S ]; then
    sed -i "$DEL_ALG_INTERFACE_S,$DEL_ALG_INTERFACE_E d" ./common/io.cpp
fi

for i in `ls`
do
    if [[ $i == alg-* ]]; then
        ALGNAME=`echo $i | awk -F- '{print $2}'`
        if [ -e ./include/$ALGNAME.h ]; then
            if [ `grep "#define ALG_REQUIRED" ./include/$ALGNAME.h | wc -l` -eq 1 ]; then
                sed -i "$ALG_OBJ_S a LIBS += $i/lib$ALGNAME.a\nLIBSFLAGS += -l$ALGNAME\n$i/lib$ALGNAME.a:\n\tcd $i; make; cd ..\n" Makefile
                sed -i "$ALG_HEADER_S a #include \"$ALGNAME.h\"" ./include/io.h
                sed -i "$ALG_INTERFACE_S a else if(!exe_$ALGNAME(seqs, alphasets, os, algo)) flag = true;" ./common/io.cpp
            fi
        fi
    fi
done