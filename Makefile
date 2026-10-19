ROOT_CONFIG ?= root-config
ROOTCLING ?= $(shell $(ROOT_CONFIG) --bindir 2>/dev/null)/rootcling
ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags 2>/dev/null)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs 2>/dev/null)
ROOT_LIBDIR := $(shell $(ROOT_CONFIG) --libdir 2>/dev/null)

MCDST_SOURCE_DIR ?= ../sthbtmaker/StRoot/McDst
MCDST_INCLUDE_DIR ?= $(MCDST_SOURCE_DIR)
MCDST_LIBRARY ?=
BUILD_DIR ?= .make
TARGET ?= converter
DICTIONARY_NAME ?= G__ConverterMcDst
CXXFLAGS ?= -O2 -g

PROJECT_CPPFLAGS = -D__ROOT__ -Iinclude -isystem "$(abspath $(MCDST_INCLUDE_DIR))"
COMPILE_FLAGS = -std=c++17 $(ROOT_CFLAGS) -Wall -Wextra -Wpedantic -MMD -MP
LINK_LIBS = $(ROOT_LIBS) -lEG -Wl,-rpath,"$(ROOT_LIBDIR)"
CORE_NAMES = OscarReader OscarConverter McDstWriter RootConverter
CORE_OBJECTS = $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(CORE_NAMES)))
MCDST_NAMES = McEvent McParticle McRun
MCDST_HEADERS = $(addprefix $(MCDST_INCLUDE_DIR)/,$(addsuffix .h,$(MCDST_NAMES)))

ifeq ($(strip $(MCDST_LIBRARY)),)
DICTIONARY_SOURCE = $(BUILD_DIR)/$(DICTIONARY_NAME).cxx
DICTIONARY_PCM = $(DICTIONARY_NAME)_rdict.pcm
MCDST_OBJECTS = $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(MCDST_NAMES))) $(BUILD_DIR)/$(DICTIONARY_NAME).o
else
MCDST_LINK = "$(abspath $(MCDST_LIBRARY))" -Wl,-rpath,"$(dir $(abspath $(MCDST_LIBRARY)))"
MCDST_DEPENDENCY = $(MCDST_LIBRARY)
endif

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.PHONY: all build test check-deps profile clean distclean help

all: $(TARGET)
build: all

$(TARGET): $(CORE_OBJECTS) $(BUILD_DIR)/converter.o $(MCDST_OBJECTS) $(MCDST_DEPENDENCY) Makefile | check-deps $(DICTIONARY_PCM)
	$(CXX) $(LDFLAGS) -o "$@" $(filter %.o,$^) $(MCDST_LINK) $(LINK_LIBS)

converter_tests: $(CORE_OBJECTS) $(BUILD_DIR)/converter_tests.o $(MCDST_OBJECTS) $(MCDST_DEPENDENCY) Makefile | check-deps $(DICTIONARY_PCM)
	$(CXX) $(LDFLAGS) -o "$@" $(filter %.o,$^) $(MCDST_LINK) $(LINK_LIBS)

test: converter_tests
	./converter_tests

$(BUILD_DIR)/%.o: src/%.cpp Makefile $(MCDST_HEADERS) | check-deps $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) -c "$<" -o "$@"

$(BUILD_DIR)/converter_tests.o: tests/converter_tests.cpp Makefile $(MCDST_HEADERS) | check-deps $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) -c "$<" -o "$@"

ifeq ($(strip $(MCDST_LIBRARY)),)
$(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(MCDST_NAMES))): $(BUILD_DIR)/%.o: $(MCDST_SOURCE_DIR)/%.cxx Makefile $(MCDST_HEADERS) | check-deps $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) -c "$<" -o "$@"

$(DICTIONARY_SOURCE): $(MCDST_HEADERS) include/ConverterMcDstLinkDef.h Makefile | check-deps $(BUILD_DIR)
	"$(ROOTCLING)" -f "$@" -D__ROOT__ -I"$(abspath $(MCDST_INCLUDE_DIR))" $(foreach header,$(MCDST_HEADERS),"$(abspath $(header))") include/ConverterMcDstLinkDef.h

$(BUILD_DIR)/$(DICTIONARY_NAME).o: $(DICTIONARY_SOURCE) Makefile | check-deps
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) -c "$<" -o "$@"

# Keep ROOT's dictionary metadata beside the executable, including when run elsewhere.
$(DICTIONARY_PCM): $(DICTIONARY_SOURCE)
	cp "$(BUILD_DIR)/$(DICTIONARY_PCM)" "$@"
endif

$(BUILD_DIR):
	mkdir -p "$@"

check-deps:
	@command -v "$(ROOT_CONFIG)" >/dev/null 2>&1 || { echo "ROOT not found: source /path/to/root/bin/thisroot.sh or set ROOT_CONFIG=/path/to/root-config." >&2; exit 1; }
ifeq ($(strip $(MCDST_LIBRARY)),)
	@test -x "$(ROOTCLING)" || { echo "rootcling not found: set ROOTCLING=/path/to/rootcling." >&2; exit 1; }
endif

profile:
	@if [ "$$(uname -s)" = Darwin ]; then echo "gprof is unavailable on macOS; use Instruments." >&2; exit 1; fi
	$(MAKE) BUILD_DIR=.make_profile TARGET=converter_profile DICTIONARY_NAME=G__ConverterMcDst_profile CXXFLAGS="$(CXXFLAGS) -pg" LDFLAGS="$(LDFLAGS) -pg" all

clean:
	rm -rf -- "$(BUILD_DIR)"
	rm -f -- "$(TARGET)" converter_tests "$(DICTIONARY_NAME)_rdict.pcm"

distclean: clean
	rm -rf -- .make_profile
	rm -f -- converter_profile G__ConverterMcDst_profile_rdict.pcm gmon.out

help:
	@echo "Run make or make test from the Converter root; CMake is optional."
	@echo "Set MCDST_SOURCE_DIR, or MCDST_INCLUDE_DIR and MCDST_LIBRARY."
	@echo "Set ROOT_CONFIG and ROOTCLING if ROOT is outside PATH."

-include $(CORE_OBJECTS:.o=.d) $(BUILD_DIR)/converter.d $(BUILD_DIR)/converter_tests.d $(MCDST_OBJECTS:.o=.d)
