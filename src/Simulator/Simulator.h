#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <math.h>

#include "../Instance/Instance.h"
#include "../Logger/Logger.h"
#include "../Round/Round.h"
#include "../SparseSet/SparseSet.h"

using namespace dataflow;

namespace simulation {

    struct State {
        int pointer;
        double usage;
        double transferRate;
    };

    struct Event {
        double time;
        int buffer;
        bool stop{false};

        bool operator<(const Event &e) const;
        bool operator==(const Event &e) const;
        bool operator<=(const Event &e) const;
        bool operator>(const Event &e) const;
        bool operator>=(const Event &e) const;
    };

    class Simulator {
    public:
        // Run with earliest stop, Rmax and Handovers are given as a target
        void runWithEarliestStopViolation(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd, double dumpRate,
                                          double startTime, double endTime, std::vector<std::vector<int>> &priority,
                                          std::vector<double> &stop, std::vector<double> &target, std::vector<double> &hTarget);
        // Run with earliest stop, Rmax and Handovers are given as a target, return first overflow
        int runWithEarliestStop(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd, double dumpRate,
                                double startTime, double endTime, std::vector<std::vector<int>> &priority,
                                std::vector<double> &stop, std::vector<double> &target, std::vector<double> &hTarget, std::vector<double> &hConstraints);
        // Run allocating the whole bandwidth to all buffers to computes lower bounds
        void runFullBandwidth(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd, double dumpRate,
                              double startTime, double endTime);

        // Run without stop
        // Rmax and handovers are not constrained
        void run(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd, double dumpRate,
                 double startTime, double endTime, std::vector<std::vector<int>> &priority);

        static const double never;

        Simulator(const Instance &instance, logger::Logger &logger);
        void initialise(std::vector<int> initialPointer, std::vector<double> initialUsage);

        void computeRates(std::vector<std::vector<int>> &priority);
        void advanceToStop(const double t);
        void setNoStop(std::vector<bool> &pNoStop);  // Set the buffer that cannot be stopped

        SparseSet<int> deactivatedBuffers;      // Indicates the stopped buffers
        std::vector<double> overflowTimes;      // Keeps track of the overflow times
        std::vector<double> capacity;           // Keeps tracks of the capacity of each buffer
        double timeOverflow;                    // Keeps track of the earliest overflow
        int bufferOverflow;                     // Keeps track of the buffer that overflows first

        double getUsage(const int b) const;
        int getPointer(const int b) const;
        double getMaxUsage(const int b) const;

    private:

        const Instance &instance;       // The instance to solve
        logger::Logger mLogger;         // Custom logger to facilitate debug
        std::vector<State> current;     // Current events during the simulation
        double now;                     // Current timepoint
        double end;                     // End time for the current downlink window
        int bufferIndex;                // Keep tracks of the buffer index of event to insert in the events
        double dumpRate;                // Current dump rate (must be updated every run)
        double precision = 0.000001;    // Precision for the values

        std::vector<double> handover;       // Handover target for the current window
        std::vector<double> handoverCt;     // Handover constraints for the current window
        std::vector<double> maxUsage;       // Max peak for each buffer
        std::vector<bool> noStop;           // Indicates if buffer can be stopped
        double getEmptyTime(const int b);
        double getEarliestStopTime(int bufferIndex, double next);
        double getFillRate(const int b);
        double getTransferRate(const int b);
        double customRound(double precision, double value);

    };
}