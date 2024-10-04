CC := g++
CFLAGS := -Wall
TARGET := hw4_1
v := 1

all: $(TARGET)

$(TARGET): f74114752_hw4_1.cpp
	$(CC) -o $@ $^ -lpthread

tsan: f74114752_hw4_1.cpp
	$(CC) -o $@ $^ -lpthread -fsanitize=thread -fPIE -pie -g

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out
