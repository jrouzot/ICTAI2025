#include "LubySequence.h"

int luby(int n) {
    int k = 1;
    // Find the largest power of 2 that is <= n + 1
    while ((1 << k) - 1 < n + 1) {
        k++;
    }
    k--;

    if ((1 << k) - 1 == n) {
        return 1 << (k - 1);
    }
    return luby(n - (1 << k) + 1);
}

LubySequence::LubySequence(int limit) {
    for (int i = 0; i < limit; i++) {
        lubySequence.push_back(luby(i + 1));
    }
}

int LubySequence::getLubyNumber(int index) {
    assert(index < lubySequence.size());
    return lubySequence[index];
}
