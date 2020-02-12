
COMPILER = g++
CXXFLAGS = -std=c++17 -O3 -DNDEBUG
MYLIBS = -I ./include
SDSLLIBS = -I ~/include -L ~/lib
SDSLFLAGS = -lsdsl -ldivsufsort -ldivsufsort64
THREADFLAGS = -pthread

SRCDIR = src
INCLUDEDIR = include
BINDIR = bin

EXECS = $(BINDIR)/fix_reference_for_indexing $(BINDIR)/construct_index $(BINDIR)/rf_mapping $(BINDIR)/contig_assembly $(BINDIR)/summary_info
all: $(EXECS)

####
$(BINDIR)/fix_reference_for_indexing:  $(SRCDIR)/fix_reference_for_indexing.cpp
	$(COMPILER) $(CXXFLAGS) $(SRCDIR)/fix_reference_for_indexing.cpp -o $(BINDIR)/fix_reference_for_indexing

$(BINDIR)/construct_index:  $(INCLUDEDIR)/genomeutil.hpp $(SRCDIR)/construct_index.cpp
	$(COMPILER) $(CXXFLAGS) $(MYLIBS) $(SDSLLIBS) $(SRCDIR)/construct_index.cpp -o $(BINDIR)/construct_index $(SDSLFLAGS)

$(BINDIR)/rf_mapping:  $(INCLUDEDIR)/genomeutil.hpp $(SRCDIR)/rf_mapping.cpp $(SRCDIR)/rf_mapping_lib.hpp
	$(COMPILER) $(CXXFLAGS) $(MYLIBS) $(SDSLLIBS) $(SRCDIR)/rf_mapping.cpp -o $(BINDIR)/rf_mapping $(SDSLFLAGS) $(THREADFLAGS)

$(BINDIR)/contig_assembly:  $(SRCDIR)/contig_assembly.cpp
	$(COMPILER) $(CXXFLAGS) $(SRCDIR)/contig_assembly.cpp -o $(BINDIR)/contig_assembly

$(BINDIR)/summary_info:  $(INCLUDEDIR)/genomeutil.hpp $(INCLUDEDIR)/bedpe_manager.hpp $(SRCDIR)/summary_info.cpp
	$(COMPILER) $(CXXFLAGS) $(MYLIBS) $(SRCDIR)/summary_info.cpp -o $(BINDIR)/summary_info

clean:
	rm -f $(EXECS)
