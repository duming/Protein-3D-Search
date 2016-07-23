#Author : Ming Du
#July 2016


CXX := g++
TITLE_NAME := CathIDX

FLAGS := -g -w -pthread

SRCS = $(wildcard *.cpp)
#OBJS := $(patsubst %.cpp,%.o,$SRCS)


main: $(SRCS)
	g++ $(FLAGS) -o $@  $^

