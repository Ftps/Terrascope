TARGET = exe
CC = g++
LIBS = -lm -lboost_iostreams -lboost_filesystem
HEAD = .
SRCS = .
INCDIRS = -I$(HEAD)
CFLAGS = -g -Wall -pedantic -O3 -fPIC -std=c++17 $(INCDIRS)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp)
MOCS = $(shell grep -l Q_OBJECT $(DEPS))
MOC_SOURCES = $(MOCS:.hpp=.moc.cpp)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp)) $(MOC_SOURCES:.cpp=.o)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f $(TARGET)
	-rm -f $(SRCS)/*.o
	-rm -f $(HEAD)/*.o
