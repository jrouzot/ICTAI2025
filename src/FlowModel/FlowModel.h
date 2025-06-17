//
// Created by Julien Rouzot on 20/05/24.
//

#pragma once

#include "../Solution/Solution.h"
#include "../Instance/Instance.h"
#include "../Logger/Logger.h"

namespace optimization {
    class FlowModel {
    public:
        FlowModel(dataflow::Instance &instance, logger::Logger &logger);

        dataflow::FillRateSolution Solve();

    private:
        dataflow::Instance mInstance;
        logger::Logger mLogger;
    };

}

