CC = mpicxx
CFLAGS = -O3 -std=c++11 -Wall
DEPS =
TARGET = MonteCarlo

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o $(TARGET)
