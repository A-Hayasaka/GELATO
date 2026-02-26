# =============================================================================
# GELATO C++ / pybind11 Build System
# =============================================================================

CXX      := g++
CXXSTD   := -std=c++17
OPTFLAGS := -O3 -march=native -DNDEBUG
WFLAGS   := -Wall -Wextra -Wno-unused-parameter
INCLUDES := -I/usr/include/eigen3
PYBIND   := $(shell python -m pybind11 --includes)
PYEXT    := $(shell python -m pybind11 --extension-suffix)

# Common compile flags for object files
CXXFLAGS := $(CXXSTD) $(OPTFLAGS) $(WFLAGS) $(INCLUDES) -fPIC $(PYBIND)

# Linker flags for shared libraries
LDFLAGS  := -shared -rdynamic

# Directories
SRCDIR   := src
OBJDIR   := build/obj
LIBDIR   := lib

# Core library sources (compiled once, shared by all modules)
CORE_SRCS := $(SRCDIR)/Air.cpp \
             $(SRCDIR)/Earth.cpp \
             $(SRCDIR)/gravity.cpp \
             $(SRCDIR)/Coordinate.cpp \
             $(SRCDIR)/iip.cpp

CORE_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CORE_SRCS))

# Pybind module sources
PYBIND_SRCS := $(SRCDIR)/pybind_USStandardAtmosphere.cpp \
               $(SRCDIR)/pybind_coordinate.cpp \
               $(SRCDIR)/pybind_dynamics.cpp \
               $(SRCDIR)/pybind_utils.cpp \
               $(SRCDIR)/pybind_IIP.cpp

PYBIND_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(PYBIND_SRCS))

# Header dependencies (for incremental builds)
WRAPPER_HDRS := $(SRCDIR)/wrapper_air.hpp \
                $(SRCDIR)/wrapper_coordinate.hpp \
                $(SRCDIR)/wrapper_utils.hpp
CORE_HDRS    := $(SRCDIR)/Air.hpp \
                $(SRCDIR)/Earth.hpp \
                $(SRCDIR)/gravity.hpp \
                $(SRCDIR)/Coordinate.hpp \
                $(SRCDIR)/iip.hpp

# Output shared libraries
TARGETS := $(LIBDIR)/USStandardAtmosphere_c$(PYEXT) \
           $(LIBDIR)/coordinate_c$(PYEXT) \
           $(LIBDIR)/dynamics_c$(PYEXT) \
           $(LIBDIR)/utils_c$(PYEXT) \
           $(LIBDIR)/IIP_c$(PYEXT)

# =============================================================================
# Build rules
# =============================================================================

all: $(TARGETS)

# Create obj directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Core object files: depend on all core + wrapper headers
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(CORE_HDRS) $(WRAPPER_HDRS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Module link rules ---
$(LIBDIR)/USStandardAtmosphere_c$(PYEXT): $(OBJDIR)/pybind_USStandardAtmosphere.o $(CORE_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LIBDIR)/coordinate_c$(PYEXT): $(OBJDIR)/pybind_coordinate.o $(CORE_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LIBDIR)/dynamics_c$(PYEXT): $(OBJDIR)/pybind_dynamics.o $(CORE_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LIBDIR)/utils_c$(PYEXT): $(OBJDIR)/pybind_utils.o $(CORE_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LIBDIR)/IIP_c$(PYEXT): $(OBJDIR)/pybind_IIP.o $(CORE_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

# =============================================================================
# Utility targets
# =============================================================================

clean:
	rm -rf $(OBJDIR)
	rm -f $(TARGETS)

rebuild: clean all

# Print build configuration (for debugging)
info:
	@echo "CXX      = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "PYEXT    = $(PYEXT)"
	@echo "CORE_OBJS= $(CORE_OBJS)"
	@echo "TARGETS  = $(TARGETS)"

.PHONY: all clean rebuild info
