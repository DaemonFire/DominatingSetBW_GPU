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
	mkdir -p bin

dothemaths: 
	./bin/main benchmarkpres/10-1.points 250 n
	./bin/main benchmarkpres/10-2.points 250 n
	./bin/main benchmarkpres/10-3.points 250 n
	./bin/main benchmarkpres/10-4.points 250 n
	./bin/main benchmarkpres/10-5.points 250 n
	./bin/main benchmarkpres/10-6.points 250 n
	./bin/main benchmarkpres/10-7.points 250 n
	./bin/main benchmarkpres/10-8.points 250 n
	./bin/main benchmarkpres/10-9.points 250 n
	./bin/main benchmarkpres/10-10.points 250 n
	./bin/main benchmarkpres/20-1.points 180 n
	./bin/main benchmarkpres/20-2.points 180 n
	./bin/main benchmarkpres/20-3.points 180 n
	./bin/main benchmarkpres/20-4.points 180 n
	./bin/main benchmarkpres/20-5.points 180 n
	./bin/main benchmarkpres/20-6.points 180 n
	./bin/main benchmarkpres/20-7.points 180 n
	./bin/main benchmarkpres/20-8.points 180 n
	./bin/main benchmarkpres/20-9.points 180 n
	./bin/main benchmarkpres/20-10.points 180 n
	./bin/main benchmarkpres/20-11-250.points 250 n
	./bin/main benchmarkpres/20-12-250.points 250 n
	./bin/main benchmarkpres/20-13-250.points 250 n
	./bin/main benchmarkpres/20-14-250.points 250 n
	./bin/main benchmarkpres/20-15-250.points 250 n
	./bin/main benchmarkpres/20-16-250.points 250 n
	./bin/main benchmarkpres/20-17-250.points 250 n
	./bin/main benchmarkpres/20-18-250.points 250 n
	./bin/main benchmarkpres/20-19-250.points 250 n
	./bin/main benchmarkpres/20-20-250.points 250 n
	./bin/main benchmarkpres/30-1.points 180 n
	./bin/main benchmarkpres/30-2.points 180 n
	./bin/main benchmarkpres/30-3.points 180 n
	./bin/main benchmarkpres/30-4.points 180 n
	./bin/main benchmarkpres/30-5.points 180 n
	./bin/main benchmarkpres/30-6.points 180 n
	./bin/main benchmarkpres/30-7.points 180 n
	./bin/main benchmarkpres/30-8.points 180 n
	./bin/main benchmarkpres/30-9.points 180 n
	./bin/main benchmarkpres/30-10.points 180 n
	./bin/main benchmarkpres/30-11-150.points 150 n
	./bin/main benchmarkpres/30-12-150.points 150 n
	./bin/main benchmarkpres/30-13-150.points 150 n
	./bin/main benchmarkpres/30-14-150.points 150 n
	./bin/main benchmarkpres/30-15-150.points 150 n
	./bin/main benchmarkpres/30-16-150.points 150 n
	./bin/main benchmarkpres/30-17-150.points 150 n
	./bin/main benchmarkpres/30-18-150.points 150 n
	./bin/main benchmarkpres/30-19-150.points 150 n
	./bin/main benchmarkpres/30-20-150.points 150 n
	./bin/main benchmarkpres/40-1.points 150 n
	./bin/main benchmarkpres/40-2.points 150 n
	./bin/main benchmarkpres/40-3.points 150 n
	./bin/main benchmarkpres/40-4.points 150 n
	./bin/main benchmarkpres/40-5.points 150 n
	./bin/main benchmarkpres/40-6.points 150 n
	./bin/main benchmarkpres/40-7.points 150 n
	./bin/main benchmarkpres/40-8.points 150 n
	./bin/main benchmarkpres/40-9.points 150 n
	./bin/main benchmarkpres/40-10.points 150 n
	./bin/main benchmarkpres/50-1.points 150 0
	./bin/main benchmarkpres/50-2.points 150 0
	./bin/main benchmarkpres/50-3.points 150 0
	./bin/main benchmarkpres/50-4.points 150 1
	./bin/main benchmarkpres/50-5.points 150 2
	./bin/main benchmarkpres/50-6.points 150 2
	./bin/main benchmarkpres/50-7.points 150 2
	./bin/main benchmarkpres/50-8.points 150 0
	./bin/main benchmarkpres/50-9.points 150 2
	./bin/main benchmarkpres/50-10.points 150 2
	./bin/main benchmarkpres/50-11-120.points 120 2
	./bin/main benchmarkpres/50-12-120.points 120 0
	./bin/main benchmarkpres/50-13-120.points 120 0
	./bin/main benchmarkpres/50-14-120.points 120 0
	./bin/main benchmarkpres/50-15-120.points 120 0
	./bin/main benchmarkpres/50-16-120.points 120 0
	./bin/main benchmarkpres/50-17-120.points 120 0
	./bin/main benchmarkpres/50-18-120.points 120 0
	./bin/main benchmarkpres/50-19-120.points 120 0
	./bin/main benchmarkpres/50-20-120.points 120 0
	./bin/main benchmarkpres/60-1.points 120 1
	./bin/main benchmarkpres/60-2.points 120 0
	./bin/main benchmarkpres/60-3.points 120 1
	./bin/main benchmarkpres/60-4.points 120 0
	./bin/main benchmarkpres/60-5.points 120 1
	./bin/main benchmarkpres/60-6.points 120 0
	./bin/main benchmarkpres/60-7.points 120 0
	./bin/main benchmarkpres/60-8.points 120 0
	./bin/main benchmarkpres/60-9.points 120 0
	./bin/main benchmarkpres/60-10.points 120 0


$(TARGET): $(OBJECTS)
	$(CC) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.cu $(HEADERS)
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

clean:
	rm -f $(TARGET)
	rm -rf obj

