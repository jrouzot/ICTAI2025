//
// Created by Julien Rouzot on 20/05/24.
//

#pragma once

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "../Logger/Logger.h"
#include "../Round/Round.h"

/*
 * This namespace is used for the instances.
 * The instances are composed of:
 *      A set of downlink windows (connection with the Deep Space Network),
 *      A set of opportunity windows (time windows with a weight for the data),
 *      A set of events (when the observations are already fixed).
 */
namespace dataflow {

    /*
     * The downlink windows have a start and end time, and a bandwidth (@rate) available.
     */
    struct Downlink {
        double start;
        double end;
        double rate;
    };

    /*
     * The opportunity windows have and start and end time and a weight (@value)
     * that indicates how interesting it is to have observations there.
     */
    struct Opportunity {
        double start;
        double end;
        double value;
    };

    /*
     * The events represent the different observations when they are fixed.
     * An event has a start time and a value that represent the new fill rate until the next event.
     * We compare events by time.
     */
    struct Event {
        double time;
        double value;

        bool operator<(const Event& other) const
        {
            return this->time < other.time;
        }

        bool operator==(const Event& other) const
        {
            return this->time == other.time && this->value == other.value;
        }
    };

    /*
     * Space missions have several instruments with different opportunity windows, events, rates and memory state.
     */
    struct Instrument {
        std::string mName;
        std::vector <Opportunity> mOpportunities;
        std::vector <Event> mEvents;
        double mMinRate;
        double mMaxRate;
        double mCapacity;
        double mInitialMemoryState;
        std::vector<double> mDataOffset;
        std::map<double, double> mDataToDumpAtTimepoint;

        // Const getter for read-only access
        const std::string& getName() const { return mName; }
        const std::vector<Opportunity>& getOpportunities() const { return mOpportunities; }
        const std::vector<Event>& getEvents() const { return mEvents; }
        void setEvent(int eventIndex, double value) { mEvents[eventIndex].value = value; }
        double getMinRate() const { return mMinRate; }
        double getMaxRate() const { return mMaxRate; }
        double getCapacity() const { return mCapacity; }
        double getInitialMemoryState() const { return mInitialMemoryState; }
        double getDataOffset(int downlinkIndex) const { return mDataOffset[downlinkIndex]; }
        const double getDataToDumpAt(double timepoint) const;
    };

    /*
     * An instance have a set of instruments and downlink windows.
     * The timepoints of an instance are the inflexion points for the transfer rates.
     * We compute the transfer rate every time something changes:
     *      When entering/exiting a downlink window,
     *      When a fill rate event occurs,
     *      When entering/exiting an opportunity window.
     */
    class Instance {
    public:
        Instance(const std::string &filename, logger::Logger &logger);

        void Print() const;

        std::vector <Instrument>::const_iterator instrumentsBegin() const;
        std::vector <Instrument>::const_iterator instrumentsEnd() const;
        std::vector <Downlink>::const_iterator downlinksBegin() const;
        std::vector <Downlink>::const_iterator downlinksEnd() const;
        std::vector <double>::const_iterator timepointsBegin() const;
        std::vector <double>::const_iterator timepointsEnd() const;

        const int &getNumInstruments() const;
        const int &getNumDownlinks() const;
        const int &getNumTimepoints() const;
        const std::string &getInstanceName() const;
        const Instrument &getInstrument(int instrumentIndex) const;
        const Downlink &getDownlink(int downlinkIndex) const;
        const std::vector <double> &getTimepoints() const;
        const double getDumpRateAt(int index) const { return mDownlinks[index].rate; };
        const double getDumpRateAtTime(double timepoint) const;
        const double getFillRateAt(int bufferIndex, double timepoint) const;
        const std::vector<int> &GetHandoverIndexes() const;
        const std::vector<int> &GetStartWindowIndexes() const;
        const std::vector<double> &GetInitialMemoryState() const;
        double getSumHandover(int windowIndex) { return mSumHandover[windowIndex]; };
        void setSumHandover(std::vector<double> &sumHandover) { mSumHandover = sumHandover; };
        double getMaxCapacity() const { return mMaxCapacity; };
        double getMemState(int t) const { return mMemState[t]; };

        // Modify the instance values to avoid big numbers that can cause the solver to crash
        void Normalize(double coef, double precision = 0.0000001); // e-6
        // Modify the instance to make only one event during the non visibility windows. This is useful to reduce the problem size.
        void Shrink();

    private:
        int mNumInstruments;                        // Number of instruments
        int mNumDownlinks;                          // Number of downlinks (we add a fake downlink at the end of the instance)
        int mNumTimePoints;                         // Number of inflexion points
        double mMaxCapacity;                        // Maximum capacity (biggest coefficient for the instance)
        std::string mInstanceName;                  // Instance name
        std::vector <Instrument> mInstruments;      // Set of instruments
        std::vector <Downlink> mDownlinks;          // Set of downlinks
        std::vector<std::vector<double>> mOffset;   // Data produced during non-visibility
        logger::Logger mLogger;                     // Logger
        std::vector<double> mTimepoints;            // The set of timepoints
        std::vector<double> mMemState;              // Global memory state at timepoint
        std::vector<int> mStartWindowIndexes;       // Indexes for the windows start
        std::vector<int> mHandoverIndexes;          // Indexes of the window ends
        std::vector<double> mInitialMemoryState;    // Initial memory state for each buffer
        std::vector<double> mSumHandover;           // The global memory state is constant when using a Round-Robin algorithm

        void InitTimepoints();                      // Init all the inflexion points
        void InitDataOffset();                      // Init all the data production during non visibility
        void InitDataToDumpAt();                    // Compute the data left to dump at each time point for each buffer
        void InitMemState();                        // Compute the global memory state at each timepoint
    };

}

