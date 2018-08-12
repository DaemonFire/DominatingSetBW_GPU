#			Stage REC Hugo Hourcade 
#		Parallel programming of algorithm solving Dominating Set based on rank-width decomposition of graphs


# Documentation: http://www.gnu.org/software/make/manual/make.html
TARGET = bin/main

CC = nvcc
CFLAGS = -g
LDFLAGS = -g

OBJECTS = $(patsubst src/%.cu, obj/%.o, $(wildcard src/*.cu))
HEADERS = $(wildcard include/*.h)

all: mkdir $(TARGET)

mkdir:
	mkdir -p obj

$(TARGET): $(OBJECTS)
	$(CC) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.cu $(HEADERS)
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

clean:
	rm -f $(TARGET)
	rm -rf obj

