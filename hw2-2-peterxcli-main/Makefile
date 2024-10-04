CC := mpicc
CFLAGS := -Wall -std=c++20 -w
TARGET := hw2_2
v := 1

all: $(TARGET)

$(TARGET):F74114752_hw2_2.c
	$(CC) -o $@ $^ $(CFLAGS)

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
	rm -f *.btr
	rm -f output/*.out

local-sync:
	mutagen sync terminate hw2-2 || exit 0
	mutagen sync create --name=hw2-2 ./ pp:/home/F74114752/hw2-2-peterxcli -m one-way-safe --ignore=hw2-2

local-jerry-sync:
	mutagen sync terminate jerry-hw2-2 || exit 0
	mutagen sync create --name=jerry-hw2-2 ./ jerry-pp:/home/F74114736/peter/hw2-2-peterxcli -m one-way-safe --ignore=hw2-2

debug: all
	mpiexec -f hosts -n 8 ./hw2_2 < input/filename/${test}.txt > tmp.out
	diff tmp.out answer/${test}.out
	rm tmp.out

remote:
	ssh  -t pp "cd hw2-2-peterxcli; bash"

remote-jerry:
	ssh  -t jerry-pp "cd /home/F74114736/peter/hw2-2-peterxcli; bash"

test: clean all
	mpiexec -f hosts -n 4 ./hw2_2 < input/filename/neib6name.txt