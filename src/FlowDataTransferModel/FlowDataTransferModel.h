#pragma once

#include "../HandoverConstraints/HandoverConstraints.h"
#include "../Instance/Instance.h"
#include "../Logger/Logger.h"
#include "../Round/Round.h"
#include "../Solution/Solution.h"

/*
 * This namespace contains the LP modeling logic.
 * We try to solve the instance and maintaining a low rmax.
 */
namespace optimization {

    enum class ObjectiveFunction {
        MIN_RMAX = 0,
        MIN_HANDOVER = 1
    };

    /*
     * We model the problem of finding the bandwidth share at each inflexion point that minimizes the rmax.
     */
    class FlowDataTransferModel {
    public:
        FlowDataTransferModel(dataflow::Instance &instance, logger::Logger &logger);
        dataflow::TransferSolution Solve(int windowStart, std::vector<double> &memoryState);

    private:
        dataflow::Instance &mInstance;
        logger::Logger mLogger;
    };

} // optimization
