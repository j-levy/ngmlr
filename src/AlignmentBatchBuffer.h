#ifndef ALIGNMENTBUFFERBATCH_H
#define ALIGNMENTBUFFERBATCH_H

#include "AlignmentBuffer.h"

class AlignmentBatchBuffer : public AlignmentBuffer 
{


	int const longReadBatchSize;

	int longReadBatchIndex;

	ReadGroup** longReadBatch;
	IntervalTree::IntervalTree<int>** readBatchCoordsTree;

public:
	AlignmentBatchBuffer(const char* const filename, const int longReadBatchSize = 20) :
            AlignmentBuffer(filename),
            longReadBatchSize(longReadBatchSize)/*, maxIntervalNumberPerKb(Config.getMaxSegmentNumberPerKb())*/
    {
         
        this->longReadBatchIndex = 0;
		// this->longReadBatch = (ReadGroup**) malloc(this->longReadBatchSize * sizeof(ReadGroup*));
		this->longReadBatch = new ReadGroup*[this->longReadBatchSize];
		this->readBatchCoordsTree = new IntervalTree::IntervalTree<int>*[this->longReadBatchSize];
    }

    ~AlignmentBatchBuffer()
    {
        delete[] this->longReadBatch;
		delete[] this->readBatchCoordsTree;
    }

    // functions related to batch processing
	void enqueueLongReadLIS(ReadGroup* group);
	void processLongReadBatchLIS(); // the read pointer is kept in memory in this->longReadBatch
    void alignSingleOrMultipleBatchIntervals(MappedRead * read, Interval const * const interval, LocationScore * tmp, Align * tmpAling, int *alignIndex, int readIndex);
    Align * alignBatchInterval(MappedRead const * const read,
		Interval const * interval, char const * const readSeq,
		size_t const readSeqLen, bool const realign, bool const fullAlignment);
    Align * computeBatchAlignment(Interval const * interval,
		int corridor, char const * const readSeq, size_t const readLength,
		int const externalQStart, int const externalQEnd, int fullReadLength,
		MappedRead const * const read, bool realign, bool const fullAlignment,
		bool const shortRead);

	int computeMappingQuality(Align const & alignment, int readLength, int readIndex);
};

#endif