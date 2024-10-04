CC := mpic++
CFLAGS := -Wall
TARGET := hw1_1
v := 1

all: $(TARGET)

$(TARGET):F74114752_hw1_1.cpp
	$(CC) -o $@ $^

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f *.btr
	rm -f output/*.out
