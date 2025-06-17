#pragma once

#include <iostream>
#include <vector>

class HandoverConstraints {

public:
    HandoverConstraints();
    // Data structure for the handover constraints added to the PL
    HandoverConstraints(const std::vector<int> &bufferSet, double target) :
            mBufferSet(bufferSet), mTarget(target) {};

    const std::vector <int> &getBufferSet() { return mBufferSet; };
    const double getTarget() { return mTarget; };

private:
    std::vector <int> mBufferSet;
    double mTarget;

};
