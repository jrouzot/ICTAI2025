#pragma once

#include "../Instance/Instance.h"
#include "../Solution/Solution.h"

namespace optimization {

    class FlowModelRelaxed {

    public:
        FlowModelRelaxed(dataflow::Instance &instance, logger::Logger &logger);
        dataflow::TransferSolution Solve(
                std::vector<std::vector<double>> &handoverConstraints,
                std::vector<std::vector<double>> &handoverFixed);

    private:
        dataflow::Instance &mInstance;
        logger::Logger mLogger;
    };
}

