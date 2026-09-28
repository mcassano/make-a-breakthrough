# Builds the CDCL search tool against an installed CaDiCaL (>= 2.0, with
# IPASIR-UP). Override CADICAL to point at another installation prefix.
CADICAL ?= $(shell brew --prefix cadical 2>/dev/null || echo /usr/local)
CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -DNDEBUG

no3_cdcl: no3_cdcl.cpp
	$(CXX) $(CXXFLAGS) -I$(CADICAL)/include -o $@ $< $(CADICAL)/lib/libcadical.a

clean:
	rm -f no3_cdcl

.PHONY: clean
