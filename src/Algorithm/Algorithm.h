#pragma once

#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "../DownlinkCount/DownlinkCount.h"
#include "../CplexFlowDataTransferModel/CplexFlowDataTransferModel.h"
#include "../HandoverConstraints/HandoverConstraints.h"
#include "../Logger/Logger.h"
#include "../LubySequence/LubySequence.h"
#include "../Options/Options.h"
#include "../Round/Round.h"
#include "../Simulator/Simulator.h"

namespace simulation {

    class Algorithm {

    public:
        Algorithm(Instance &i, logger::Logger &logger);

        double playDownlinkViolation(const int i, std::vector<std::vector<int>> &priority,
                                     std::vector<double> &handoverTargets, std::vector<double> &handoverConstraints);
        void playDownlinkEarliestStop(const int i, std::vector<std::vector<int>> &priority,
                                      std::vector<double> &handoverTargets, std::vector<double> &handoverConstraints);
        void playDownlinkWithoutStop(const int i, std::vector<std::vector<int>> &priority);
        void playDownlinkFullBandwidth(const int i);
        void playSolution(std::vector<std::vector<std::vector<int>>> &priority);
        void playSolutionFullBandwidth();
        void computeHandoverUb();


        double solve(double timeLimit);
        bool checkSolution(double rmax);


        const std::vector<std::vector<int>>& getBestPriority() const { return bestPriority; }

        // Performances indicators
        int nbLoops = 0;
        int nbSingleWindow = 0;
        int nbSimulation = 0;
        int nbCallsMove = 0;
        int nbFailsMove = 0;
        double timeTotal = 0;
        double timePL = 0.0;
        double timeMove = 0.0;

    private:
        Instance &instance;
        logger::Logger &mLogger;
        Simulator simulator;
        std::vector<int> buffers;
        std::vector<std::vector<double>> stop;
        std::vector<std::vector<double>> handover;
        std::vector<std::vector<double>> handoverLb;
        std::vector<std::vector<double>> handoverUb;
        std::vector<std::vector<double>> maxUsage;
        std::vector<std::vector<double>> usage;
        std::vector<std::vector<int>> pointer;
        std::vector<double> sumHandover;
        std::vector<Event> events;
        std::vector<std::vector<Event>::const_iterator> window;
        std::vector<double> dump;
        std::vector<double> target;
        std::vector<double> swTarget;
        double rmaxUpperBound;
        double rmaxLowerBound;
        double rmax;
        double singleWindowMinViolation;
        double precision = 0.00001;

        // Dummy priority assignment for each window
        std::vector<std::vector<std::vector<int>>> dummyPriority;
        // Dummy priority assignment for one window
        std::vector<std::vector<int>> dummyPriorityWindow;
        // Dummy handover for each window
        std::vector<std::vector<double>> dummyHandovers;
        // Best priority assignment for a single window
        std::vector<std::vector<int>> bestPriority;
        // Best priority assignment in a different format
        std::vector<int> flatPriority;
        // Overflowing buffers for minViolation
        std::vector<int> overflowingBuffersMinViolation;
        // Overflowing buffers
        std::vector<int> overflowingBuffers;
        // Priority assigment of best solution
        std::vector<std::vector<std::vector<int>>> prioritySolution;
        // Temporary priority assignment
        std::vector<std::vector<std::vector<int>>> priorityTmp;
        // Handovers of best solution
        std::vector<std::vector<double>> handoverSolution;
        // Current buffer to prevent stop
        std::vector<bool> noStop;

        void setTargetRmax(double rmax);        // Update the target to capacity * rmax
        void setSwTargetRmax(double rmax);      // Update the swTarget to capacity * rmax
        bool postHandoverConstraint(int windowIndex);
        void runDownlinkCount(std::vector<double> &initialMemoryState, int nbRuns);
        bool relaxHandovers(int windowIndex);
        void setNoStop(int windowIndex);
        bool singleWindow(int windowIndex, double targetRmax, std::vector<double> &currentHandoverTargets, std::vector<double> &currentHandoverConstraints);
        void initHandoverConstraints();
        void updateHandoverTargets();
    };
}