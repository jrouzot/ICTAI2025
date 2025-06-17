#pragma once

#include <cassert>
#include <stdlib.h>
#include <vector>

class LubySequence {
public:
    LubySequence(int limit);
    int getLubyNumber(int index);
private:
    std::vector<int> lubySequence;
};
