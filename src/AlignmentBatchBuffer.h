#ifndef ALIGNMENTBUFFERBATCH_H
#define ALIGNMENTBUFFERBATCH_H

#include "AlignmentBuffer.h"

class AlignmentBatchBuffer : public AlignmentBuffer 
{


	int const longReadBatchSize;

	int longReadBatchIndex;

	ReadGroup** longReadBatch;

public:
	AlignmentBatchBuffer(const char* const filename, const int longReadBatchSize = 1000) :
            AlignmentBuffer(filename),
            longReadBatchSize(longReadBatchSize)/*, maxIntervalNumberPerKb(Config.getMaxSegmentNumberPerKb())*/
    {
         
        this->longReadBatchIndex = 0;
		// this->longReadBatch = (ReadGroup**) malloc(this->longReadBatchSize * sizeof(ReadGroup*));
		this->longReadBatch = new ReadGroup*[this->longReadBatchSize];

    }

    ~AlignmentBatchBuffer()
    {
        delete this->longReadBatch;
    }

    // functions related to batch processing
	void enqueueLongReadLIS(ReadGroup* group);
	void processLongReadBatchLIS(); // the read pointer is kept in memory in this->longReadBatch
    void alignSingleOrMultipleBatchIntervals(MappedRead * read, Interval const * const interval, LocationScore * tmp, Align * tmpAling, int & alignIndex);

};

#endif