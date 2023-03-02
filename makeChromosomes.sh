#!/bin/bash

# this is not a makefile - just command line compile

#g++ -ggdb -Wall -std=c++11 -o pchrom pchromosomes.cpp
g++ -O3 -Wall -std=c++11 -o pchrom parameters.cpp random.cpp pchromosomes.cpp

