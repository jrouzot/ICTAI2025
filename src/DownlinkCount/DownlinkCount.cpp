#include <iostream>
#include <random>

#include "DownlinkCount.h"

DownlinkCount::DownlinkCount(dataflow::Instance &instance, double rmaxTarget) :
                            mInstance(instance), mTarget(rmaxTarget) {}

std::vector<std::vector<int>> DownlinkCount::run(int windowIndex, std::vector<double> &initialMemoryState) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mTarget, mTarget/4);
    double target = dist(gen);
    if(target < 0.1) {
        target = 0.1;
    }

    mCurrentUsage.resize(mInstance.getNumInstruments(), 0.0);
    std::vector<bool> overflow(mInstance.getNumInstruments(), false);
    mPriority = {};
    for(int b{0}; b < mInstance.getNumInstruments(); ++b) {
        mPriority.push_back({});
    }
    int currentPriorityAssignment = mInstance.getNumInstruments()-1;
    bool windowOverflow;
    int cnt = 0;

    for(int b{0}; b < mInstance.getNumInstruments(); ++b) {
        mCurrentUsage[b] = initialMemoryState[b];
    }

    for(int w{windowIndex}; w < mInstance.getNumDownlinks(); ++w) {
        windowOverflow = false;
        for(int b{0}; b < mInstance.getNumInstruments(); ++b) {
            if(not(overflow[b])) {
                mCurrentUsage[b] += mInstance.getInstrument(b).getDataOffset(windowIndex);
                mCurrentUsage[b] += mInstance.getInstrument(b).getDataToDumpAt(mInstance.getDownlink(windowIndex).start);
                if(mCurrentUsage[b] >= mInstance.getInstrument(b).getCapacity() * target) {
                    windowOverflow = true;
                    overflow[b] = true;
                    mPriority[currentPriorityAssignment].emplace_back(b);
                    cnt++;
                }
            }
        }
        if(cnt == mInstance.getNumInstruments()) {
            break;
        }
        if(windowOverflow) {
            --currentPriorityAssignment;
            assert(currentPriorityAssignment >= 0);
        }
    }
    // Add the remaining buffers in the worst priority pool
    for(int b{0}; b < mInstance.getNumInstruments(); ++b) {
        if (not(overflow[b])) {
            mPriority[currentPriorityAssignment].emplace_back(b);
        }
    }
    // Remove empty priorities
    while (!mPriority.empty() && mPriority.front().empty()) {
        mPriority.erase(mPriority.begin());
    }

    return mPriority;
}



