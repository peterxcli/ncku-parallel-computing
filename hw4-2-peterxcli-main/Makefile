CC := g++
CFLAGS := -Wall
TARGET := hw4_2
v := 1

all: $(TARGET)

$(TARGET):f74114752_hw4_2.cpp
	$(CC) -o $@ $^ -lpthread -std=c++11

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out
