CC := g++
CFLAGS := -Wall
TARGET := final2
v := 1

all: $(TARGET)

$(TARGET):f74114752_final2.cpp
	$(CC) -o $@ $^ -fopenmp

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
