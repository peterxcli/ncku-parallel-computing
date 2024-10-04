CC := mpicxx
CFLAGS := -Wall
TARGET := hw6
v := 1

$(TARGET):f74114752_hw6.cpp
	$(CC) -o $@ $^ -lpthread -std=c++11

local-sync:
	mutagen sync terminate hw6 || exit 0
	mutagen sync create --name=hw6 ./ pp:/home/F74114752/hw6-peterxcli --ignore=hw6

clean:
	rm -f $(TARGET)
	rm -f output/*.out

run:
	mpiexec -n 8 -f hosts ./$(TARGET) < ${input} > ${output}
