#include "AlignmentBatchBuffer.h"

void AlignmentBatchBuffer::enqueueLongReadLIS(ReadGroup* group)
{
	this->longReadBatch[this->longReadBatchIndex] = group;
	this->longReadBatchIndex++;
	if (this->longReadBatchIndex >= this->longReadBatchSize)
	{
		this->processLongReadBatchLIS();
		this->longReadBatchIndex = 0;
	}
}

void AlignmentBatchBuffer::processLongReadBatchLIS()
{
	
}