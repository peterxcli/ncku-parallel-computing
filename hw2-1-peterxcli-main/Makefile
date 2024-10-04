CC := mpic++
CFLAGS := -Wall
TARGET := hw2_1
v := 1

all: $(TARGET)

$(TARGET):F74114752_hw2_1.cpp
	$(CC) -o $@ $^

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f output/*.out
	rm -f *.btr


local-sync:
	mutagen sync terminate hw2-1 || exit 0
	mutagen sync create --name=hw2-1 ./ pp:/home/F74114752/hw2-1-peterxcli

debug: all
	mpiexec -f hosts -n 8 ./hw2_1 < input/filename/${test}.txt > tmp.out
	diff tmp.out answer/${test}.out
	rm tmp.out
