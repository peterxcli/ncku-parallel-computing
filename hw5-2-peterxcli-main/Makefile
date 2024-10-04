CC := g++
CFLAGS := -Wall
TARGET := hw5_2
v := 1

all: $(TARGET)

$(TARGET):f74114752_hw5_2.cpp
	$(CC) -o $@ $^ -Xpreprocessor -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -std=c++17 || $(CC) -o $@ $^ -fopenmp -std=c++17

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out

local-sync:
	mutagen sync terminate hw5-2 || exit 0
	mutagen sync create --name=hw5-2 ./ pp:/home/F74114752/hw5-2-peterxcli --ignore=hw5-2

