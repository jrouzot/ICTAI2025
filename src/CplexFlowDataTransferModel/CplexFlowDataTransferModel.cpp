#include "CplexFlowDataTransferModel.h"

CplexFlowDataTransferModel::CplexFlowDataTransferModel(dataflow::Instance &instance, logger::Logger &logger)
: mInstance(instance), mLogger(logger), env(), model(env) {}

dataflow::TransferSolution CplexFlowDataTransferModel::Solve() {

    /***********************
     *                     *
     * -    VARIABLES    - *
     *                     *
     ***********************/

    /*
     * The transfer rate variables keeps tracks of the new transfer rate at each inflexion point, for each buffer
     */
    std::vector<std::vector<IloNumVar>> transferRatesVars({});
    transferRatesVars.resize(mInstance.getNumTimepoints());
    for (auto t(0); t < mInstance.getNumTimepoints(); ++t) {
        // Only add decision variable for inflexion points inside visibility windows
        if (mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]) > 0 or t == mInstance.getNumTimepoints() - 2) {
            double timeLag = mInstance.getTimepoints()[t + 1] - mInstance.getTimepoints()[t];
            transferRatesVars[t].resize(mInstance.getNumInstruments());
            for (auto b(0); b < mInstance.getNumInstruments(); ++b) {
                IloNumVar transferVar(env, 0.0, mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]) * timeLag);
                transferRatesVars[t][b] = transferVar;
            }
        } else {
            transferRatesVars[t] = {};
        }
    }

    /*
     * The memory variables keep track of the memory state at each inflexion point for each buffer
     */
    std::vector<std::vector<IloNumVar>> memoryVars({});
    memoryVars.resize(mInstance.getNumTimepoints());
    for (int t(0); t < mInstance.getNumTimepoints(); ++t) {
        // Only add decision variable for inflexion points inside visibility windows
        if ((mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]) > 0 or t == mInstance.getNumTimepoints() - 2)
            or std::find(mInstance.GetHandoverIndexes().begin(), mInstance.GetHandoverIndexes().end(), t) !=
               mInstance.GetHandoverIndexes().end()) {
            memoryVars[t].resize(mInstance.getNumInstruments());
            for (auto b(0); b < mInstance.getNumInstruments(); ++b) {
                IloNumVar memoryVar(env, 0.0, mInstance.getInstrument(b).getCapacity());
                memoryVars[t][b] = memoryVar;
            }
        } else {
            memoryVars[t] = {};
        }
    }

    /*
     * Auxiliary variable for the rmax
     */
    IloNumVar alpha(env, 0);
/*    IloNumVar alpha(env, 0.0, 1.0);*/


    /*************************
     *                       *
     * -     OBJECTIVE     - *
     *                       *
     *************************/

    /*IloObjective obj = IloMinimize(env, mInstance.getMaxCapacity() * alpha - 10e-3 * sumTotalTransfer);*/
/*    IloObjective obj = IloMinimize(env, 1 * alpha - 10-2 * sumTotalTransfer);*/
    IloObjective obj = IloMinimize(env, alpha);
    model.add(obj);


    /*************************
     *                       *
     * -    CONSTRAINTS    - *
     *                       *
     *************************/

    /*
     * Memory start
     */
    for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
        double memoryStart =
                mInstance.getInstrument(b).getInitialMemoryState() + mInstance.getInstrument(b).getDataOffset(0);
        model.add(memoryVars[0][b] == memoryStart);
    }

    /*
     * Bandwidth share
     */
    for (auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
        if (transferRatesVars[t].size()) {
            double timeLag = mInstance.getTimepoints()[t + 1] - mInstance.getTimepoints()[t];
            double bandwidth = mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]);
            // sum(transferVar at t)
            IloExpr sumTransfer(env);
            for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                sumTransfer += transferRatesVars[t][b];
            }
            model.add(sumTransfer <= bandwidth * timeLag);
        }
    }

    /*
     * Memory consistency: At each time point the memory sum is constant
     */
    for (auto t(0); t < mInstance.getNumTimepoints(); ++t) {
        assert(memoryVars.size());
        IloExpr sumMem(env);
        for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            sumMem += memoryVars[t][b];
        }
        model.add(sumMem == mInstance.getMemState(t));
    }

    /*
     * Memory consistency (same window)
     */
    int windowIndex = 0;
    for (auto t{1}; t < mInstance.getNumTimepoints(); ++t) {
        if (transferRatesVars[t - 1].size()) {
            double timeLag = mInstance.getTimepoints()[t] - mInstance.getTimepoints()[t - 1];
            assert(timeLag > 0.0);
            for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                double fillRate = mInstance.getFillRateAt(b, mInstance.getTimepoints()[t - 1]);
                // mem(i,j) = mem(i,j-1) + tr(i-1,j) * timeLag
                model.add(memoryVars[t - 1][b] - transferRatesVars[t - 1][b] + fillRate * timeLag == memoryVars[t][b]);
            }
        } else {
            ++windowIndex;
            for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                double sum = mInstance.getInstrument(b).getDataOffset(windowIndex);
                model.add(memoryVars[t-1][b] + sum == memoryVars[t][b]);
            }
        }
    }

/*    *//*
     * Memory consistency (non visibility)
     *//*
    for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
        int prevIndex;
        int windowIndex;
        for (int w{1}; w < mInstance.GetStartWindowIndexes().size(); ++w) {
            prevIndex = mInstance.GetHandoverIndexes()[w - 1];
            windowIndex = mInstance.GetStartWindowIndexes()[w];
            double sum = mInstance.getInstrument(b).getDataOffset(w);
            model.add(memoryVars[prevIndex][b] + sum == memoryVars[windowIndex][b]);
        }
    }*/

    /*
     * Handover consistency (the sum of handovers is constant so all the bandwidth is used)
     */
/*    int currentWindowIndex = 0;
    for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
        // Add the constraints only for the handover timepoints
        if(mInstance.GetHandoverIndexes()[currentWindowIndex] == t) {
            IloExpr sumMem(env);
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                sumMem += memoryVars[t][b];
            }
            model.add(sumMem == mInstance.getSumHandover(currentWindowIndex));
            ++currentWindowIndex;
        }
    }*/

    /*
     *  mem / capacity <= alpha
     */
    for (auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
        for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            model.add(1 * memoryVars[t][b] / mInstance.getInstrument(b).getCapacity() <= alpha);
/*                model.add(memoryVars[t][b] / mInstance.getInstrument(b).getCapacity()  <= alpha);*/
        }
    }

    // Create the model
    cplex = IloCplex(model);
    cplex.setParam(IloCplex::Param::Emphasis::Numerical, IloTrue);

    // Désactiver les logs
    cplex.setOut(env.getNullStream()); // Désactiver la sortie standard
    cplex.setWarning(env.getNullStream()); // Désactiver les avertissements

    // Save current time
    double startTime = cplex.getCplexTime();

    // Solve the model
    if (cplex.solve()) {

        double solveTime = cplex.getCplexTime() - startTime;

        mLogger.Debug("Optimal solution found:");
        mLogger.Debug("Objective value = ", cplex.getValue(alpha));
        mLogger.Debug("Solve time = ", solveTime);

        // Fetch handover values
        std::vector<std::vector<double>> handover;
        handover.resize(mInstance.getNumDownlinks());
        int currentWindowIndex = 0;
        for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {

            if (mInstance.GetHandoverIndexes()[currentWindowIndex] == t) {
                handover[currentWindowIndex].resize(mInstance.getNumInstruments());
                double sum = 0.0;
                for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                    handover[currentWindowIndex][b] = cplex.getValue(memoryVars[t][b]);
                    sum += handover[currentWindowIndex][b];
                }
                assert(sum - 10e-3 <= mInstance.getMemState(t) and mInstance.getMemState(t) <= sum + 10e-3);
                ++currentWindowIndex;
            }
        }

/*
        // Print variables
        std::cout << "Variables:\n";
        for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                std::cout << "m" << t << b << "=";
                std::cout << cplex.getValue(memoryVars[t][b]) << "\n";
            }
        }
*/

        dataflow::TransferSolution solution(cplex.getValue(alpha), handover, solveTime);
        return solution;

    } else {
        mLogger.Error("No solution found in the PL !");
        exit(1);
    }
}
