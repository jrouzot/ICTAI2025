#pragma once

#include "../Instance/Instance.h"
#include "../Simulator/Simulator.h"

class DownlinkCount {
public:
    DownlinkCount(Instance &instance, double rmaxTarget);
    std::vector<std::vector<int>> run(int windowIndex, std::vector<double> &initialMemoryState);

private:
    Instance mInstance;
    double mTarget;
    std::vector<double> mCurrentUsage;
    std::vector<std::vector<int>> mPriority;
};
