CC := mpicxx
CFLAGS := -Wall
TARGET := mid4
v := 1

all: $(TARGET)

$(TARGET):f74114752_mid4.cpp
	$(CC) -o $@ $^

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out
