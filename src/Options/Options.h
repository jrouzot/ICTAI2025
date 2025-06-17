#pragma once

enum class BufferSelectionHeuristic {
    OVERFLOW = 0,
    MOST_MEM = 1,
    RANDOM = 2
};

enum class ValueSelectionHeuristic {
    MAX_VIOLATION = 0,
    MIN_VIOLATION = 1,
    RANDOM = 2
};

/*
 * Represents the options that can be given by the users from the command line.
 */
struct Options {
    BufferSelectionHeuristic bufferSelectionHeuristic;
    ValueSelectionHeuristic valueSelectionHeuristic;
};