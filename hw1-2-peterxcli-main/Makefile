CC := mpic++
CFLAGS := -Wall
TARGET := hw1_2
v := 1

all: $(TARGET)

$(TARGET):F74114752_hw1_2.cpp
	$(CC) -o $@ $^ -w

judge: clean all 
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out
	rm -f *btr
