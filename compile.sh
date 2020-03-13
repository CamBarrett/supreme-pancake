
C_COMPILER="g++"
C_INCLUDES="-I/usr/include"
C_FLAGS="-O2 -W -g"


if [ -e myutil.o ]; then
    rm myutil.o
fi
${C_COMPILER} ${C_FLAGS} ${C_INCLUDES} -c myutil.C

if [ -e my_main_sequence.o ]; then
    rm my_main_sequence.o
fi
${C_COMPILER} ${C_FLAGS} ${C_INCLUDES} -c my_main_sequence.C

if [ -e DWD_pop_BP_ugriz ]; then
    rm DWD_pop_BP_ugriz
fi
${C_COMPILER} ${C_FLAGS} ${C_INCLUDES} -o DWD_pop_BP_ugriz DWD_pop_BP_ugriz.C myutil.o my_main_sequence.o


