bin_PROGRAMS = CreateBloom FindExons

CreateBloom_CPPFLAGS = -I$(top_srcdir)/Common

CreateBloom_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/StringUtil.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
        Common/nthist.hpp \
        Common/nthash.hpp \
        Common/ntHashIterator.hpp \
        Common/BloomFilter.hpp \
        CreateBloom.cpp

FindExons_CPPFLAGS = -I$(top_srcdir)/Common

FindExons_SOURCES = \
        Common/FastaReader.cpp \
        Common/FastaReader.h \
        Common/FastaConcat.h \
        Common/IOUtil.h \
        Common/StringUtil.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/Sequence.cpp \
        Common/Sequence.h \
        Common/Options.cpp \
        Common/Options.h \
        Common/pstream.h \
        Common/nthist.hpp \
        Common/nthash.hpp \
        Common/ntHashIterator.hpp \
        Common/BloomFilter.hpp \
        FindExons.cpp

dist_bin_SCRIPTS= MakeSplicegraph.py

dist_doc_DATA = \
	ChangeLog \
	COPYRIGHT \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh
