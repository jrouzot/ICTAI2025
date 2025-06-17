#pragma once

#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <iostream>

#include "../Instance/Instance.h"
#include "../Logger/Logger.h"
#include "../Solution/Solution.h"

class CplexFlowDataTransferModel {
public:
    CplexFlowDataTransferModel(dataflow::Instance &instance, logger::Logger &logger);
    ~CplexFlowDataTransferModel() { env.end(); };
    dataflow::TransferSolution Solve();

private:
    dataflow::Instance &mInstance;
    logger::Logger mLogger;

    IloEnv env;
    IloModel model;
    IloCplex cplex;
};
