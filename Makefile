#Author : Ming Du
#July 2016


CXX := g++
TITLE_NAME := CathIDX

FLAGS := -g -w -pthread
INCLUDE_PATH := -Ilib


#$(patsubst %.cpp,$(OBJ_DIR)/%.o,$

LIB_DIR := $(sort $(dir $(wildcard lib/*/*)))
OBJ_DIR := obj

MAIN_OBJ := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(wildcard *.cpp))
LIB_OBJ :=  $(wildcard $(patsubst  %,%*.cpp,$(LIB_DIR)))



#default: 
#	@echo $(LIB_OBJ)

executable : $(MAIN_OBJ) 
	$(CXX) $(flags) $(INCLUDE_PATH) -o $(TITLE_NAME) $^


#$(LIB_OBJ) : 
#$(OBJ_DIR)/lib/%.o : lib/%.cpp
#	$(CXX) $(flags) $(INCLUDE_PATH)  -c $< -o $@

$(MAIN_OBJ) :
$(OBJ_DIR)/%.o : %.cpp
	$(CXX) $(flags)  $(INCLUDE_PATH)  -c $< -o $@

