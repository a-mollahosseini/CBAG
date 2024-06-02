# Boost library
BOOST = Path_To_Your_Boost_Here # Change to your Boost library

# C++ compiler
CXX = g++-11 # Change to your C++ compiler
CXXFLAGS = -lstdc++ -std=c++17 -O2

SOURCES = $(wildcard Code/*.cpp Code/BetweennessApprox/*.cpp)

CBAG: $(SOURCES)
	$(CXX) -I $(BOOST) $(SOURCES) $(CXXFLAGS) -o CBAG

clean:
	rm CBAG
