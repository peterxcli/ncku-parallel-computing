CC := g++
CFLAGS := -Wall
TARGET := hw5_1
v := 1

all: $(TARGET)

$(TARGET):f74114752_hw5_1.cpp
	$(CC) -o $@ $^ -fopenmp

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out

local-sync:
	mutagen sync terminate hw5-1 || exit 0
	mutagen sync create --name=hw5-1 ./ pp:/home/F74114752/hw5-1-peterxcli --ignore=hw5-1