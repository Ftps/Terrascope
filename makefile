TARGET = exe
CC = g++
MOC = moc-qt5
LIBS = -Lqcustomplot -lm -lqcustomplot -lboost_system -lboost_filesystem -lpthread -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -lQt5Widgets -lQt5Gui -lQt5Core -lQt5PrintSupport
HEAD = ./Headers
SRCS = ./Source
INCDIRS = -I/usr/include/qt -I/usr/include/qt/QtWidgets -I/usr/include/qt/QtCore -I/usr/include/qt/QtGui -I$(HEAD)
CFLAGS = -Wall -pedantic -O3 -fPIC -std=c++17 $(INCDIRS)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp)
MOCS = $(shell grep -l Q_OBJECT $(DEPS))
MOC_SOURCES = $(MOCS:.hpp=.moc.cpp)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp)) $(MOC_SOURCES:.cpp=.o)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

%.moc.cpp: %.hpp
	$(MOC) $(INCDIRS) $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f $(TARGET)
	-rm -f $(SRCS)/*.o
	-rm -f $(HEAD)/*.o
