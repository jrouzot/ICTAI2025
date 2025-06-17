//
// Created by Julien Rouzot on 20/05/24.
//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Instance.h"
#include "../Simulator/Simulator.h"

namespace dataflow {

    /*
     * Parse instance file to an Instance object.
     */
    Instance::Instance(const std::string &filename, logger::Logger &logger) :
            mLogger(logger) {

        mLogger.Debug("Initializing Instance...");

        mInstanceName = filename;

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open instance file: " + filename);
        }

        int numInstruments;
        std::string flag;

        file >> numInstruments;
        mNumInstruments = numInstruments;
        file >> flag;

        mInstruments.resize(numInstruments);

        assert(flag == "instruments");

        int current;
        std::vector <std::string> names(numInstruments);
        std::vector <double> minRates(numInstruments);
        std::vector <double> maxRates(numInstruments);
        std::vector <double> capacities(numInstruments);
        std::vector <double> memoryStart(numInstruments);

        for(int i(0); i < numInstruments; ++i) {
            file >> names[i] >> minRates[i] >> maxRates[i] >> memoryStart[i] >>  capacities[i];
        }

        for(int i(0); i < numInstruments; ++i) {
            memoryStart[i] = memoryStart[i];
            capacities[i] = capacities[i];
        }

        mInitialMemoryState = memoryStart;
        int numDownlinks;

        file >> numDownlinks;
        mNumDownlinks = numDownlinks+1;     // We add a dummy window at end
        file >> flag;

        assert(flag == "downlinks");

        mDownlinks.resize(mNumDownlinks);

        double start;
        double end;
        double value;

        for(auto i(0); i < numDownlinks; ++i) {
            file >> current;
            assert(current == i);
            file >> start >> end >> value;
            mDownlinks[i] = {start, end, value};
        }

        std::string tmp;
        std::string name;
        double startOpportunity;
        double endOpportunity;
        double weightOpportunity;
        std::vector<std::vector<Opportunity>> opportunities({});

        for(int i(0); i < numInstruments; ++i) {
            int numOpportunities;
            std::vector<Opportunity> opportunity;
            file >> numOpportunities;
            file >> tmp >> tmp >> name;

            assert(name == names[i]);

            opportunity.resize(numOpportunities);

            for (auto j(0); j < numOpportunities; ++j) {
                file >> current;
                assert(current == j);
                file >> startOpportunity;
                file >> endOpportunity;
                file >> weightOpportunity;
                opportunity[j] = {startOpportunity, endOpportunity, weightOpportunity};
            }
            opportunities.push_back(opportunity);
        }

        std::vector<std::vector<Event>> eventList({});
        double t;
        double v;

        for(int i(0); i < numInstruments; ++i) {
            int numEvents;
            std::vector<Event> events;
            file >> numEvents;
            file >> tmp >> tmp >> name;

            events.resize(numEvents);

            for(auto j(0); j < numEvents; ++j) {
                file >> t;
                file >> v;
                events[j] = {
                        t,
                        v
                };
            }
            eventList.push_back(events);
        }

        for(int i(0); i < numInstruments; ++i) {
            mInstruments[i] = {
                    names[i],
                    opportunities[i],
                    eventList[i],
                    minRates[i],
                    maxRates[i],
                    capacities[i],
                    memoryStart[i]
            };
        }

        // Add a dummy window to end the instance
        double maxTime{0.0};
        for(int i{0}; i < numInstruments; ++i) {
            if(eventList[i].size() and eventList[i][eventList[i].size()-1].time >= maxTime) {
                maxTime = eventList[i][eventList[i].size()-1].time;
            }
        }
        mDownlinks[mNumDownlinks-1] = {maxTime+1, maxTime+2, 0};

        // Find max capacity
        mMaxCapacity = 0.0;
        for(int b{0}; b < numInstruments; ++b) {
            if(mInstruments[b].getCapacity() > mMaxCapacity) {
                mMaxCapacity = mInstruments[b].getCapacity();
            }
        }

        mLogger.Debug("Successfully fetched instance");
        mTimepoints = {};
        mHandoverIndexes = {};
        this->Normalize(10e-6, 10e-9);
        this->InitDataOffset();
        this->Shrink();
        this->InitTimepoints();
        mNumTimePoints = mTimepoints.size();
        this->InitMemState();
        // Compute the sum of the fill until end of downlink at each timepoint
        this->InitDataToDumpAt();
/*        this->Print();*/
    }

    /*
     * Print the Instance object in a nice way.
     */
    void Instance::Print() const {
        std::cout << "[INSTANCE]: " << mInstanceName << std::endl;
        std::cout << "\n";
        for(int i(0); i < mDownlinks.size(); ++i) {
            std::cout << "\tDownlink " << i+1 << ": " << mDownlinks[i].start << " " << mDownlinks[i].end << " " << mDownlinks[i].rate << std::endl;
        }
        std::cout << "\n";
        for(const Instrument& instrument : mInstruments) {
            std::cout << "\t" << instrument.mName << ": " << std::endl;
            std::cout << "\t\t" << "min rate: " << instrument.mMinRate << std::endl;
            std::cout << "\t\t" << "max rate: " << instrument.mMaxRate << std::endl;
            std::cout << "\t\t" << "capacity: " << instrument.mCapacity << std::endl;
            std::cout << "\t\t" << "initial memory state: " << instrument.mInitialMemoryState << std::endl;
            std::cout << "\t\t" << "opportunities: " << std::endl;
            for(auto opportunity : instrument.mOpportunities) {
                std::cout << "\t\t\t" << opportunity.start << " " << opportunity.end << " " << opportunity.value << std::endl;
            }
            std::cout << "\n";
            std::cout << "\t\t" << "events: " << std::endl;
            for(auto event : instrument.mEvents) {
                std::cout << "\t\t\t" << event.time << " " << event.value << std::endl;
            }
        }
        std::cout << "Timepoints: " << std::endl;
        std::cout << "\t";
        for(auto t : mTimepoints) {
            std::cout << t << " ";
        }
        std::cout << std::endl;

        std::cout << "Handover timepoints indexes: " << std::endl;
        std::cout << "\t";
        for(auto t : mHandoverIndexes) {
            std::cout << t << " ";
        }
        std::cout << std::endl;

        for(int b{0}; b < mNumInstruments; ++b) {
            std::cout << "DataOffset for buffer " << b << ": " << std::endl;
            std::cout << "\t";
            for (auto offset: mInstruments[b].mDataOffset) {
                std::cout << offset << " ";
            }
            std::cout << std::endl;
        }

        for(auto b{0}; b < mNumInstruments; ++b) {
            std::cout << "dataToDumpAt for buffer: " << b << std::endl;
            std::cout << "\t";
            for(auto t : mTimepoints) {
                std::cout << mInstruments[b].getDataToDumpAt(t) << " ";
            }
            std::cout << std::endl;
        }
    }

    /*
     * Get the number of instruments
     */
    const int& Instance::getNumInstruments() const {
        return mNumInstruments;
    }

    /*
     * Get the number of instruments
     */
    const int& Instance::getNumDownlinks() const {
        return mNumDownlinks;
    }

    /*
     * Get the number of timepoints
     */
    const int& Instance::getNumTimepoints() const {
        return mNumTimePoints;
    }

    /*
     * Get the timepoint vector
     */
    const std::vector <double>& Instance::getTimepoints() const {
        return mTimepoints;
    }

    /*
     * Get the handover indexes vector
     */
    const std::vector <int>& Instance::GetHandoverIndexes() const {
        return mHandoverIndexes;
    }

    /*
     * Get the window start indexes vector
     */
    const std::vector <int>& Instance::GetStartWindowIndexes() const {
        return mStartWindowIndexes;
    }

    /*
     * Get the initial memory state
     */
    const std::vector <double>& Instance::GetInitialMemoryState() const {
        return mInitialMemoryState;
    }

    /*
     * Get the bandwidth available at a given timepoint
     */
    const double Instance::getDumpRateAtTime(double timepoint) const {
        for(auto d : mDownlinks) {
            if(timepoint < d.start) {
                return 0.0;
            }
            if(timepoint < d.end) {
                return d.rate;
            }
        }
        return 0.0;
    }

    /*
     * Get the fill rate for a given buffer at a given timepoint
     */
    const double Instance::getFillRateAt(int bufferIndex, double timepoint) const {
        double rate = 0.0;
        for(auto e : mInstruments[bufferIndex].getEvents()) {
            if(timepoint < e.time) {
                return rate;
            }
            rate = e.value;
        }
        return rate;
    }

    /*
     * Get the instruments begin iterator.
     */
    std::vector <Instrument>::const_iterator Instance::instrumentsBegin() const {
        return mInstruments.cbegin();
    }

    /*
     * Get the instruments end iterator.
     */
    std::vector <Instrument>::const_iterator Instance::instrumentsEnd() const {
        return mInstruments.cend();
    }

    /*
     * Get the downlinks begin iterator.
     */
    std::vector <Downlink>::const_iterator Instance::downlinksBegin() const {
        return mDownlinks.cbegin();
    }

    /*
     * Get the downlinks end iterator.
     */
    std::vector <Downlink>::const_iterator Instance::downlinksEnd() const {
        return mDownlinks.cend();
    }

    /*
     * Get the timepoints start iterator
     */
    std::vector <double>::const_iterator Instance::timepointsBegin() const {
        return mTimepoints.cbegin();
    }

    /*
     * Get the timepoints end iterator
     */
    std::vector <double>::const_iterator Instance::timepointsEnd() const {
        return mTimepoints.cend();
    }

    /*
     * Get a downlink by index.
     */
    const Downlink& Instance::getDownlink(int downlinkIndex) const {
        assert(downlinkIndex < mDownlinks.size());
        return mDownlinks[downlinkIndex];
    }

    /*
     * Get instrument by index.
     */
    const Instrument& Instance::getInstrument(int instrumentIndex) const {
        assert(instrumentIndex < mInstruments.size());
        return mInstruments[instrumentIndex];
    }

    /*
     * Get instance name.
     */
    const std::string& Instance::getInstanceName() const {
        return mInstanceName;
    }

    /*
     * Get data to dump at timepoint
     */
    const double Instrument::getDataToDumpAt(double timepoint) const {
        auto it = mDataToDumpAtTimepoint.lower_bound(timepoint);
        // Found the element
        if (it != mDataToDumpAtTimepoint.end() && it->first == timepoint) {
            return it->second;
        }

        // No valid previous element
        if (it == mDataToDumpAtTimepoint.begin()) {
            std::cerr << "No key less than or equal to " << timepoint << std::endl;
            exit(1);
        }

        // Compute data to dump
        double timeAfter = it->first;
        double valueAfter = it->second;
        --it;
        double timeBefore = it->first;
        double valueBefore = it->second;
        double rate = (valueBefore - valueAfter) / (timeAfter - timeBefore);

        return valueBefore - rate * (timepoint - timeBefore);
    }


    /*
     * Set a new timepoint for each event during downlink window @donwlinkIndex
     */
    void Instance::InitTimepoints() {

        for(const auto & instrument : mInstruments) {
            for(const auto & event : instrument.getEvents()) {
                // Insert the time of the fill rate event
                mTimepoints.push_back(event.time);
            }
        }

        for(auto & downlink : mDownlinks) {
            mTimepoints.emplace_back(downlink.start);
            mTimepoints.emplace_back(downlink.end);
        }

        // Sort the timepoint by date
        std::sort(mTimepoints.begin(), mTimepoints.end());
        // Remove duplicates
        mTimepoints.erase(std::unique(mTimepoints.begin(), mTimepoints.end()), mTimepoints.end());
        // Save number of timepoints
        mNumTimePoints = mTimepoints.size();

        // Save the downlink end timepoints indexes
        int currentDownlinkIndex = 0;
        for(int t{0}; t < mNumTimePoints; ++t) {
            if(mTimepoints[t] == mDownlinks[currentDownlinkIndex].start) {
                mStartWindowIndexes.emplace_back(t);
            }
            if(mTimepoints[t] == mDownlinks[currentDownlinkIndex].end) {
                mHandoverIndexes.emplace_back(t);
                ++currentDownlinkIndex;
            }
        }
    }

    /*
     * Compute the global memory state at each time point
     */
    void Instance::InitMemState() {
        int currentWindow = 0;
        mMemState.resize(mTimepoints.size());
        double sum = 0.0;
        for(int b{0}; b < mNumInstruments; ++b) {
            sum += getInstrument(b).getInitialMemoryState() + getInstrument(b).getDataOffset(currentWindow);
        }
        mMemState[0] = sum;
        for(int t{0}; t < mNumTimePoints - 1; ++t) {
            int deltaT = mTimepoints[t+1] - mTimepoints[t];
            int time = mTimepoints[t];
            if(getDumpRateAtTime(mTimepoints[t]) > 0) {
                for(int b{0}; b < mNumInstruments; ++b) {
                    double fillRate = getFillRateAt(b, mTimepoints[t]);
                    sum += fillRate * deltaT;
                }
            } else {
                ++currentWindow;
                for(int b{0}; b < mNumInstruments; ++b) {
                    sum += getInstrument(b).getDataOffset(currentWindow);
                }
            }
            sum -= getDumpRateAtTime(mTimepoints[t]) * deltaT;
            sum = std::max(0.0, sum);
            mMemState[t+1] = sum;
        }
    }


    /*
     * Compute the data produced during non visibility and save it as a constant for each buffer window
     */
    void Instance::InitDataOffset() {
        for(int b{0}; b < mNumInstruments; ++b) {
            mInstruments[b].mDataOffset.resize(mNumDownlinks);
            for(int w{0}; w < mNumDownlinks; ++w) {
                double sum = 0.0;
                double currentValue = 0.0;
                double startTime = 0.0;
                if(w > 0) {
                    startTime = mDownlinks[w-1].end;
                }
                double currentTime = startTime;

                // Loop through all buffer events
                for(auto e : mInstruments[b].getEvents()) {

                    // We only check events before the current downlink window
                    if (e.time < mDownlinks[w].start) {
                        // We only add value if the event is inside the non visibility
                        if(e.time >= startTime) {
                            sum += (e.time - currentTime) * currentValue;
                            currentTime = e.time;
                        }
                        currentValue = e.value;
                        // We add the last segment
                    } else {
                        sum += (mDownlinks[w].start - currentTime) * currentValue;
                        break;
                    }
                }

                // Update data offset
                mInstruments[b].mDataOffset[w] = sum;
            }
        }
    }

    /*
     * Compute the data left to dump at each time point for each buffer
     */
    void Instance::InitDataToDumpAt() {
        for (auto b{0}; b < getNumInstruments(); ++b) {
            int windowIndex = 0;
            for(auto t{timepointsBegin()}; t < timepointsEnd(); ++t) {
                // Find the downlink window of current timepoint
                while(getDownlink(windowIndex).end < *t) {
                    ++windowIndex;
                }
                double currentEndWindowTime = getDownlink(windowIndex).end;
                auto currentEvent = mInstruments[b].getEvents().begin();
                auto timepointIterator = t;
                while((currentEvent+1) != mInstruments[b].getEvents().end() and (currentEvent+1)->time <= *timepointIterator) {
                    ++currentEvent;
                }
                double sum = 0.0;
                while((currentEvent+1) <= mInstruments[b].getEvents().end()) {
                    auto time = std::max(currentEvent->time, *timepointIterator);
                    if((currentEvent+1) == mInstruments[b].getEvents().end() or (currentEvent+1)-> time >= currentEndWindowTime) {
                        sum += (currentEndWindowTime - time) * currentEvent->value;
                        break;
                    }
                    sum += ((currentEvent+1)->time - time) * currentEvent->value;
                    ++currentEvent;
                    ++timepointIterator;
                }
                mInstruments[b].mDataToDumpAtTimepoint[*t] = sum;
            }
        }
    }

    /*
     * Multiply the instance value by @coef, and round to match @precision
     */
    void Instance::Normalize(double coef, double precision) {
        for(auto i{0}; i < mNumInstruments; ++i) {
            mInstruments[i].mCapacity = std::round(mInstruments[i].mCapacity * coef / precision) * precision;
            mInstruments[i].mInitialMemoryState = std::round(mInstruments[i].mInitialMemoryState * coef / precision) * precision;
        }
        for(auto w{0}; w < mNumDownlinks; ++w) {
            mDownlinks[w].rate = std::round(mDownlinks[w].rate * coef / precision) * precision;
        }
        for(auto i{0}; i < mNumInstruments; ++i) {
            for(auto j{0}; j < mInstruments[i].getEvents().size(); ++j) {
                mInstruments[i].setEvent(j ,std::round(mInstruments[i].getEvents()[j].value * coef / precision) * precision);
            }
        }
    }


    /*
     * Eliminates all non visibility buffer events and add a start end end event for each window
     */
    void Instance::Shrink() {
        for(int i{0}; i < mNumInstruments; ++i) {
            std::vector<Event> events = {};
            bool beforeStart = true;
            int windowIndex = 0;
            double nextWindowTime = mDownlinks[windowIndex].start;
            double currentFill = 0.0;
            // Add an end event
            mInstruments[i].mEvents.push_back({mDownlinks[mNumDownlinks- 1].start, 0});
            auto e = mInstruments[i].getEvents().begin();
            while(e < mInstruments[i].getEvents().end()) {

                if (beforeStart and e->time < nextWindowTime)  {
                    /*mLogger.Debug("Before start, ignoring event");*/
                    currentFill = e->value;
                    ++e;
                } else if(not(beforeStart) and e->time < nextWindowTime) {
                    /*mLogger.Debug("Before end, adding event");*/
                    events.push_back({e->time, e->value});
                    currentFill = e->value;
                    ++e;
                } else if(beforeStart and e->time == nextWindowTime) {
                    /*mLogger.Debug("At start, going to downlink end");*/
                    nextWindowTime = mDownlinks[windowIndex].end;
                    beforeStart = false;
                    currentFill = e->value;
                } else if(beforeStart and e->time > nextWindowTime) {
                    /*mLogger.Debug("After start, adding start event, going to downlink end");*/
                    events.push_back({nextWindowTime, currentFill});
                    nextWindowTime = mDownlinks[windowIndex].end;
                    beforeStart = false;
                    /*currentFill = e->value;*/
                } else if (not(beforeStart) and e->time >= nextWindowTime) {
                    /*mLogger.Debug("After end, adding end event, going to next downlink start");*/
                    ++windowIndex;
                    nextWindowTime = mDownlinks[windowIndex].start;
                    beforeStart = true;
                    /*currentFill = e->value;*/
                }
            }
            while(windowIndex < mNumDownlinks - 1) {
                ++windowIndex;
                events.push_back({mDownlinks[windowIndex].start, 0});
            }
            mInstruments[i].mEvents = events;
        }
    }
}
