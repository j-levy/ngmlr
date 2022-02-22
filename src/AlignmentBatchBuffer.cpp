#include "AlignmentBatchBuffer.h"

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <float.h> //Eclipse
#include <limits.h>

#include "Timing.h"
#include "SequenceProvider.h"
#include "StrippedSW.h"
#include "OutputReadBuffer.h"

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


// J.L. change function signature to admit a batch of all of this 
void AlignmentBatchBuffer::alignSingleOrMultipleBatchIntervals(MappedRead * read, Interval const * const interval, LocationScore * tmp, Align * tmpAling, int & alignIndex) {

	int readSeqLen = interval->onReadStop - interval->onReadStart;
	auto readPartSeq = extractReadSeq(readSeqLen, interval, read);

	if (readPartSeq != 0) {
		Align * align = alignInterval(read, interval, readPartSeq.get(), readSeqLen, false, false);
		if (align != 0) {
			if (align->Score > 0.0f) {
				int svType = SV_NONE;

				if (Config.getSmallInversionDetection() || Config.getLowQualitySplit()) {
					Interval * leftOfInv = new Interval();
					Interval * rightOfInv = new Interval();
					svType = SV_UNKNOWN;
					bool inversionAligned = false;
					svType = detectMisalignment(align, interval, readPartSeq.get(), leftOfInv, rightOfInv, read);

					if (svType != SV_NONE) {
						int mq = computeMappingQuality(*align, read->length);
						int assumedSvType = svType;
						svType = realign(svType, interval, leftOfInv, rightOfInv, read, tmpAling, alignIndex, tmp, mq);
					} else {
						verbose(0, true, "No SV detected!");
					}
					if (leftOfInv != 0) {
						delete leftOfInv;
						leftOfInv = 0;
					}
					if (rightOfInv != 0) {
						delete rightOfInv;
						rightOfInv = 0;
					}
				}

				if (svType == SV_NONE) {
					verbose(0, true, "No SV was detected. Keeping single alignment!");
					/**********************************************/
					// No inversion detected
					/**********************************************/
					if (satisfiesConstraints(align, read->length)) {
						align->MQ = computeMappingQuality(*align, read->length);
						align->clearNmPerPosition();

						tmpAling[alignIndex] = *align;
						delete align;
						align = 0;

						tmp[alignIndex].Location.m_Location = interval->onRefStart + tmpAling[alignIndex].PositionOffset;					//- (corridor >> 1); handled in computeAlingment
						tmp[alignIndex].Location.setReverse(interval->isReverse);
						tmp[alignIndex].Score.f = tmpAling[alignIndex].Score;

						tmpAling[alignIndex].mappedInterval = getIntervalFromAlign(&tmpAling[alignIndex], &tmp[alignIndex], alignIndex, read->length);

						read->Calculated += 1;
						alignIndex += 1;
					} else {
						align->clearBuffer();
						align->clearNmPerPosition();
						delete align;
						align = 0;
						verbose(0, true, "Alignment did not satisfiesConstraints");
					}
				} //else {
					if (align != 0) {
						align->clearBuffer();
						align->clearNmPerPosition();
						delete align;
						align = 0;
					}
//				}
			} else {
				if (align != 0) {
					align->clearBuffer();
					align->clearNmPerPosition();
					delete align;
					align = 0;
				}

				verbose(0, true, "Alignment failed");
			}
		}
	} else {
		Log.Message("Extracting read sequence failed for read %s. Please report this on https://github.com/philres/ngmlr", read->name);
	}
}



void AlignmentBatchBuffer::processLongReadBatchLIS() {

    fprintf(stderr, "\tprocess batch LETZGO\n");
    
    // initialize arrays for batch processing.
    int *nIntervalsBatch = new int[this->longReadBatchSize];


    int i = 0;
    for(i = 0; i < this->longReadBatchIndex; i++)
    {
        ReadGroup * group = this->longReadBatch[i];  
        
        MappedRead * read = group->fullRead;

        int maxAnchorNumber = 10000;
        int const maxNumScores = 1000;

        /**
         * Get all mapping positions from anchors (non overlapping 256bp parts of reads)
         * and convert them to Anchor objects
         * TODO: remove fixed length of 100
         */
        int maxCandidatesPerReadPart = 100;

        Timer overallTmr;
        overallTmr.ST();

        if (read->length <= readPartLength) {
            Log.Message("Should not be here. Short read in long read pipeline.");
            WriteRead(group->fullRead, false);
        }

        /*
        * Parts of read that were aligned plus MQ of subalignments
        */
        std::vector<IntervalTree::Interval<int> > treeIntervals;

    //	Log.Message("Processing LongReadLIS: %d - %s (length %d)", read->ReadId, read->name, read->length);

        verbose(0, true, "Processing LongReadLIS: %d - %s (length %d)", read->ReadId, read->name, read->length);
        verbose(0, true, "Anchor list:");

        Anchor * anchorsFwd = new Anchor[maxAnchorNumber];
        int anchorFwdIndex = 0;


    //	float * bestScores = new float[group->readNumber];
    //	int bestScoresIndex = 0;
    //	for (int j = 0; j < group->readNumber; ++j) {
    //		MappedRead * part = group->reads[j];
    //		if(part->numScores() > 0) {
    //			bestScores[bestScoresIndex++] = part->Scores[0].Score.f;
    //		}
    //	}
    //
    //	std::sort(bestScores, bestScores + bestScoresIndex);
    //
    //	for(int i = 0; i < bestScoresIndex; ++i) {
    //		verbose(1, false, "%f, ", bestScores[i]);
    //	}
    //	verbose(0, true, "");
    //	float const minScore = bestScores[(int)(bestScoresIndex * 0.8f)];
    //	verbose(1, true, "Min score: %f", minScore);
    //
    //	delete[] bestScores;
    //	bestScores = 0;

    //	Log.Message("Read: %s", group->fullRead->name);
    //	Log.Message("Parts: %d", group->readNumber);
    //	for (int j = 0; j < group->readNumber; ++j) {
    //		MappedRead * part = group->reads[j];
    //
    //		int positionOnRead = j * readPartLength;
    //
    //		Log.Message("\tPart %d:", j);
    //		Log.Message("\t\tName: %s", part->name);
    //		Log.Message("\t\tMappings: %d", part->numScores());
    //		Log.Message("\t\tQuality: %d", part->mappingQlty);
    //
    //		for (int k = 0; k < part->numScores(); ++k) {
    //			Log.Message("\t\t\t%d, %llu, %f", positionOnRead, part->Scores[k].Location.m_Location, part->Scores[k].Score.f );
    //		}
    //
    //	}


        /**
         * Collect sub-read mapping locations and
         * build interval tree with read locations
         * plus mapping quality
         */
        for (int j = 0; j < group->readNumber; ++j) {
            MappedRead * part = group->reads[j];

            int positionOnRead = j * readPartLength;

            verbose(1, false, "%d: ", positionOnRead);

            if (part->numScores() < maxNumScores) {
                if (part->numScores() > 0) {

                    /**
                     * Add mapped read parts + mapping quality.
                     * Used to estimate MQ for read(part)s
                     */
                    treeIntervals.push_back(IntervalTree::Interval<int>(positionOnRead, positionOnRead + readPartLength, part->mappingQlty));

                    for (int k = 0; k < part->numScores(); ++k) {

                        if (stdoutPrintScores) {
                            printf("%f\n", part->Scores[k].Score.f);
                        }

                        /**
                         * Anchor is valid and will be used
                         */
                        if (anchorFwdIndex >= (maxAnchorNumber - 1)) {
                            verbose(0, true, "Anchor array too small - reallocating.");
                            maxAnchorNumber = maxAnchorNumber * 2;
                            Anchor * anchorsTmp = new Anchor[maxAnchorNumber];
                            memcpy(anchorsTmp, anchorsFwd, anchorFwdIndex * sizeof(Anchor));
                            delete[] anchorsFwd;
                            anchorsFwd = anchorsTmp;
                        }
                        Anchor & anchor = anchorsFwd[anchorFwdIndex++];

                        anchor.score = part->Scores[k].Score.f;
                        anchor.isReverse = part->Scores[k].Location.isReverse();
                        anchor.type = DP_STATUS_OK;
                        anchor.isUnique = part->numScores() == 1; // || anchor.score > minScore;

                        /**
                         * It would be best to convert reads or read parts that map to the negative strand
                         * Immediately to plus strand.
                         * If somebody tells me how to do this properly, I'll happily change this.
                         * For now Anchors can be on the plus or minus strand!
                         * Problem: Transforming coordinates from negative anchors to plus strand
                         * is not possible without knowing if the read originates from the plus
                         * or the minus strand. This is difficult for reads with few matching anchors
                         * and reads that originate from e.g. an inverted translocation.
                         */
                        anchor.onRead = positionOnRead;
                        anchor.onRef = part->Scores[k].Location.m_Location;
                        if (anchor.isReverse) {
                            printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, anchor.onRead, anchor.onRead + readPartLength, part->Scores[k].Location.m_Location + readPartLength, part->Scores[k].Location.m_Location, part->Scores[k].Score.f,
                                    part->Scores[k].Location.isReverse(),
                                    DP_TYPE_UNFILTERED, anchor.isUnique ? DP_STATUS_OK : DP_STATUS_LOWSCORE);
                        } else {
                            printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, anchor.onRead, anchor.onRead + readPartLength, part->Scores[k].Location.m_Location, part->Scores[k].Location.m_Location + readPartLength, part->Scores[k].Score.f,
                                    part->Scores[k].Location.isReverse(),
                                    DP_TYPE_UNFILTERED, anchor.isUnique ? DP_STATUS_OK : DP_STATUS_LOWSCORE);
                        }

                        if (k < 3) {
                            double const K = 1 / 3;
                            double  const lambda = 1.098612;
                            double const n = 256;
                            double const m = 3000000000;
                            double const Y = part->Scores[k].Score.f;
                            double pVal = 1 - exp(-K * n * m * exp(-lambda * Y));
                            verbose(0, false, "%f (%f) at %llu, ", part->Scores[k].Score.f, pVal, part->Scores[k].Location.m_Location);
                        } else if (k == 3) {
                            verbose(0, false, "... (%d)", part->numScores());
                        }
                    }
                    verbose(0, true, "");
                } else {
                    verbose(0, true, "no hits found");
                    printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, positionOnRead, positionOnRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
                }
            } else {
                verbose(0, true, "too many hits found");
                printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, positionOnRead, positionOnRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
            }
        }

        Anchor * anchorsRev = new Anchor[maxAnchorNumber];
        int anchorRevIndex = 0;

        /**
         * Tree contains all mapped read parts + MQ
         */
        readCoordsTree = new IntervalTree::IntervalTree<int>(treeIntervals);

        /**
         * Build HSP from sub read mappings (cLIS)
         * - Sort sub-reads by read position
         * - Compute max. number of HSP (based on read length)
         * - Run cLIS algorithm on reference positions to retrieve HSP
         * - If isUnique, build Interval from sub-read set and compute regression
         */
        int nIntervals = 0;
        Interval * * intervals = getIntervalsFromAnchors(nIntervals, anchorsFwd, anchorFwdIndex, anchorsRev, anchorRevIndex, group->fullRead);

        verbose(0, true, "\nIntervals after cLIS:");
        for (int i = 0; i < nIntervals; ++i) {
            verbose(0, "", intervals[i]);
        }
        verbose(0, true, "");

        std::sort(intervals, intervals + nIntervals, sortIntervalsInSegment);

        std::vector<IntervalTree::Interval<Interval *> > intervalList;

        verbose(0, true, "\n\nBuilding segments:\n");
        /**
         * A segment is a list of intervals that are located in on alignment corridor
         * Intervals that are contained in others will be deleted. All the others
         * will be added to a segment
         */
        int const maxMappedSegementCount = nIntervals + 1;
        MappedSegment * segments = new MappedSegment[maxMappedSegementCount];
        size_t segementsIndex = 0;
        for (int i = 0; i < nIntervals; ++i) {
            Interval * interval = intervals[i];
            bool intervalProcessed = false;
            verbose(0, "Current interval: ", interval);

            for (int j = 0; j < segementsIndex && !intervalProcessed; ++j) {
                verbose(1, true, "Checking segment %d", j);

                for (int k = 0; k < segments[j].length && !intervalProcessed; ++k) {
                    Interval * processedInterval = segments[j].list[k];

                    verbose(2, "Comparing to interval: ", processedInterval);
                    if (isContained(interval, processedInterval)) {
                        // Interval is contained in already processed interval and therefore ignored
                        intervalProcessed = true;
                        verbose(3, true, "Is contained. Deleting interval.");
                        delete intervals[i];
                        intervals[i] = 0;
                        interval = 0;
                    } else {
                        if (isCompatible(interval, processedInterval)) {
                            verbose(3, true, "Is compatible");
                            // Interval fits corridor of segment
                            if (segments[j].length < segments[j].maxLength) {
                                segments[j].list[segments[j].length++] = interval;
                                intervalProcessed = true;
                                intervalList.push_back(IntervalTree::Interval<Interval *>(interval->onReadStart, interval->onReadStop, interval));
                            }
                        } else {
                            verbose(3, true, "Not contained and not compatible");
                        }
                    }
                }

            }

            if (!intervalProcessed) {
                // Creating new segment for interval
                verbose(1, true, "Creating new segment %d", segementsIndex);
                if (segementsIndex < maxMappedSegementCount) {
                    segments[segementsIndex].list[0] = interval;
                    segments[segementsIndex++].length = 1;
                    intervalList.push_back(IntervalTree::Interval<Interval *>(interval->onReadStart, interval->onReadStop, interval));
                } else {
                    Log.Error("Could not create new segment (%d, %d) for read %s", segementsIndex, maxMappedSegementCount, read->name);
                }
            }

        }

        /**
         * All intervals are either deleted or added to segments.
         */
        delete[] intervals;
        intervals = 0;

        // Join segments
        if (pacbioDebug) {
            Log.Message("\nSegments identified: %d", segementsIndex);
            for(int i = 0; i < segementsIndex; ++i) {
                Log.Message("Segment %d contains %d intervals", i, segments[i].length);
            }
        }

        IntervalTree::IntervalTree<Interval *> * intervalsTree = new IntervalTree::IntervalTree<Interval *>(intervalList);

        verbose(0, true, "\n\nMerging segments:\n");
        // Final interval list
        intervals = new Interval *[maxSegmentCount + 1];
        nIntervals = 0;

        Interval * * delIntervals = new Interval *[maxSegmentCount + 1];
        int nDelIntervals = 0;


        /**
         * Joins segments to full length alignments
         *  - try to identify all segments that fall into an alignment corridor
         *  - merge segments that are separated by regions with high error rate
         *  - merge segments that are separated by small indels (only if both segments are long enough to compensate the score penalty of the deletion)
         *  - Extend unmerged segments to cover the full read an compensate for missed subread mappings
         *  - Don't merge segments if separated by a read part that maps to a different place on the genome or is inverted (balanced translocation)
         */
        verbose(0, true, "Joining segments to intervals:");
        for (int i = 0; i < segementsIndex; ++i) {

            // Sort intervals by position on read
            std::sort(segments[i].list, segments[i].list + segments[i].length, sortIntervalsInSegment);

            if (pacbioDebug) {
                Log.Message("Segment %d:", i);
                for(int j = 0; j < segments[i].length; ++j) {
                    verbose(0, "Interval: ", segments[i].list[j]);
                }
            }

            Interval * lastInterval = segments[i].list[0];
            extendIntervalStart(lastInterval, 2 * readPartLength);
            bool isFirstInterval = true;
            Interval * currentInterval = 0;

            for (int j = 1; j < segments[i].length; ++j) {
                currentInterval = segments[i].list[j];
                verbose(1, "Last: ", lastInterval);
                verbose(1, "Current: ", currentInterval);
                if (isSameDirection(currentInterval, lastInterval)) {
                    verbose(2, true, "Same direction.");
                    loc dupLength = 0;
                    if (!isDuplication(currentInterval, lastInterval, dupLength)) {

                        if(gapOverlapsWithInterval(lastInterval, currentInterval, intervalsTree, read)) {
                            /***
                             * Possible translocation
                             */
                            verbose(2, true, "Overlap found in other corridor.");
                            verbose(2, true, "Saving last interval. Using current to go on.");

                            if (isFirstInterval) {
                                extendToReadStart(lastInterval, read->length, intervalsTree, read);
                                isFirstInterval = false;
                            }
                            //TODO: close gap on read (imporve)
                            extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
                            extendIntervalStart(currentInterval, 2 * readPartLength);
                            //Add lastInterval
                            intervals[nIntervals++] = lastInterval;
                            lastInterval = currentInterval;
                        } else {
                            /**
                             * Insertion, deletion or gap
                             */
                            REAL const corridorSize = std::min(4096, std::min(currentInterval->lengthOnRead(), lastInterval->lengthOnRead()));
                            verbose(2, true, "IsContained in corridor of %f.", corridorSize);
                            if (canSpanDeletionInsertion(currentInterval, lastInterval, corridorSize) && !spansChromosomeBorder(currentInterval, lastInterval)) {
                                /**
                                 *  Deletion or insertion small enough for alignment without split
                                 */
                                verbose(2, true, "Not a duplication, not overlapping with other corridor. Merging current and last interval. Deleting current interval.");
                                lastInterval = mergeIntervals(lastInterval, currentInterval);
                                segments[i].list[j]->isProcessed = true;
                                delIntervals[nDelIntervals++] = segments[i].list[j];
                            } else {
                                /**
                                 * Deletion, insertion or gap
                                 */
                                verbose(2, true, "Diagonal offset between intervals too big (max %f) or spans chromosome border.", corridorSize);
                                verbose(2, true, "Saving last interval. Using current to go on.");

                                if (isFirstInterval) {
                                    extendToReadStart(lastInterval, read->length, intervalsTree, read);
                                    isFirstInterval = false;
                                }
                                closeGapOnRead(lastInterval, currentInterval, read->length);
                                extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
                                extendIntervalStart(currentInterval, 2 * readPartLength);
                                //Add lastInterval
                                intervals[nIntervals++] = lastInterval;
                                lastInterval = currentInterval;
                            }
                        }
                    } else {
                        /**
                         * Duplication
                         */
                        verbose(2, true, "Covers possible duplication.");
                        verbose(2, true, "Saving last interval. Using current to go on.");
                        if (isFirstInterval) {
                            extendToReadStart(lastInterval, read->length, intervalsTree, read);
                            isFirstInterval = false;
                        }
                        closeGapOnRead(lastInterval, currentInterval, read->length);
                        // Extension necessary since border can be wrong even if there is no gap on the read
                        int maxExtend = std::min(std::max(currentInterval->onReadStart - lastInterval->onReadStop + (int)dupLength, 0), 2 * readPartLength);
                        verbose(2, true, "Duplication stats: %d - %d + %d, %d => %d", currentInterval->onReadStart, lastInterval->onReadStop, (int)dupLength, 2 * readPartLength, maxExtend);
    //					int maxExtend = 2 * readPartLength;
                        extendIntervalStop(lastInterval, maxExtend, read->length);
                        extendIntervalStart(currentInterval, maxExtend);

                        //Add lastInterval
                        intervals[nIntervals++] = lastInterval;
                        lastInterval = currentInterval;
                    }
                } else {
                    /**
                     * Inversion
                     */
                    verbose(2, true, "Not in same direction.");
                    verbose(2, true, "Saving last interval. Using current to go on.");
                    if (isFirstInterval) {
                        extendToReadStart(lastInterval, read->length, intervalsTree, read);
                        isFirstInterval = false;
                    }
                    closeGapOnRead(lastInterval, currentInterval, read->length);
                    extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
                    extendIntervalStart(currentInterval, 2 * readPartLength);
                    //Add lastInterval
                    intervals[nIntervals++] = lastInterval;
                    lastInterval = currentInterval;
                }
            }
            if (isFirstInterval) {
                extendToReadStart(lastInterval, read->length, intervalsTree, read);
                isFirstInterval = false;
            }
            extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
            verbose(2, true, "Extending last interval to read end", intervalsTree);
            extendToReadStop(lastInterval, read->length, intervalsTree, read);
            verbose(1, "Last: ", lastInterval);
            verbose(2, true, "Saving last interval.");
            intervals[nIntervals++] = lastInterval;
        }

        delete[] segments;
        segments = 0;

        delete intervalsTree;
        intervalsTree = 0;

        for (int i = 0; i < nDelIntervals; ++i) {
            delete delIntervals[i];
            delIntervals[i] = 0;
        }
        delete[] delIntervals;
        delIntervals = 0;
        nDelIntervals = 0;

        verbose(0, true, "Joined intervals:");
        for (int i = 0; i < nIntervals; ++i) {
            verbose(1, "", intervals[i]);
        }
        verbose(0, true, "Sorting intervals by read start");
        std::sort(intervals, intervals + nIntervals, sortIntervalsInSegment);


        if (nIntervals > 0) {
            Interval * lastInterval = intervals[0];
            for (int i = 1; i < nIntervals; ++i) {
                Interval * currentIntervl = intervals[i];

                if (currentIntervl->anchorLength > 1) {
                    verbose(1, "a:", lastInterval);
                    verbose(1, "b: ", currentIntervl);
                    if (!isCompatible(lastInterval, currentIntervl) && getDistanceOnRead(lastInterval, currentIntervl) > 0 && (currentIntervl->anchorLength > 2 || lastInterval->anchorLength > 2)) {
                        verbose(1, true, "Closing gap between:");
                        closeGapOnRead(lastInterval, currentIntervl, read->length);
                    } else {
                        verbose(1, true, "Skip");
                    }
                }

                if (currentIntervl->anchorLength > 1 || lastInterval->anchorLength == 1) {
                    lastInterval = currentIntervl;
                }
            }
        }

        /**
         * Sort intervals by score
         * Important because, we trust intervals with higher
         * score more. Thus we align them first. Aligned intervals
         * are considered fixed. Therefore, all unangliend intervals will
         * be trimmed in order to not overlap with fixed intervals.
         */
        verbose(0, true, "Sorting intervals by score");
        std::sort(intervals, intervals + nIntervals, sortIntervalsByScore);

        int readBpCovered = 0;
        verbose(0, true, "\nFinal intervals:");
        for (int i = 0; i < nIntervals; ++i) {
            verbose(1, "", intervals[i]);
            printDotPlotLine(read->ReadId, read->name, intervals[i]->onReadStart, intervals[i]->onReadStop, intervals[i]->onRefStart, intervals[i]->onRefStop, intervals[i]->score, intervals[i]->isReverse,
            DP_TYPE_SEQMENTS_CONS + i, DP_STATUS_OK);
            readBpCovered += intervals[i]->lengthOnRead();
        }

        float aligned = readBpCovered * 1.0f / read->length;
        verbose(0, true, "Intervals cover %.2f%% of read", aligned * 100.0f);
        bool mapped = (Config.getMinResidues() < 1.0) ? aligned > Config.getMinResidues() : readBpCovered > Config.getMinResidues();
        if (!mapped) {
            verbose(0, true, "Clearing intervals -> read unmapped");
            for (int i = 0; i < nIntervals; ++i) {
                if (intervals[i] != 0) {
                    delete intervals[i];
                    intervals[i] = 0;
                }
            }
            delete[] intervals;
            intervals = 0;
            nIntervals = 0;

        }
    // }


        /**
         * Align final intervals to the reference
         */
        if (nIntervals != 0) {

            // Since we don't know how many segments of the read we have to align in the
            // end we need temp arrays to store them
            // TODO: remove fixed upper limit
            LocationScore * tmpLocationScores = new LocationScore[nIntervals * 4];
            Align * tmpAlingments = new Align[nIntervals * 4];
            int nTempAlignments = 0;

            read->Calculated = 0;

            Timer tmr;
            if (pacbioDebug) {
                tmr.ST();
            }
            /**
             * Aligning intervals
             */
            for (int i = 0; i < nIntervals; ++i) {
                Interval * currentInterval = intervals[i];

                /**
                 * Adjust unanglined intervals to not overlap with already
                 * aligned intervals
                 */
                verbose(0, "Aligning interval: ", currentInterval);
                for (int j = 0; j < nTempAlignments; ++j) {
                    Interval * alignedInterval = tmpAlingments[j].mappedInterval;
                    verbose(1, "Check overlap with: ", alignedInterval);

                    int overlap = getOverlapOnRead(currentInterval, alignedInterval);
                    verbose(1, true, "Overlap: %d", overlap);
                    if (overlap > 0 && overlap < currentInterval->lengthOnRead() * 0.95f) {
                        verbose(0, true, "Adjusting interval");

                        if (currentInterval->onReadStart < alignedInterval->onReadStart) {
                            shortenIntervalEnd(currentInterval, overlap);
                        } else {
                            shortenIntervalStart(currentInterval, overlap);
                        }

                    }
                }
                verbose(0, "New interval: ", currentInterval);

                /**
                 * Convert intervals on reverse strand to forward strand
                 */
                if (currentInterval->onRefStart > currentInterval->onRefStop) {
                    loc tmp = currentInterval->onRefStart;
                    currentInterval->onRefStart = currentInterval->onRefStop;
                    currentInterval->onRefStop = tmp;
                }

                verbose(0, "Aligning interval: ", currentInterval);
                if (!Config.getSkipalign()) {
                    // J.L. change this
                    alignSingleOrMultipleIntervals(read, currentInterval, tmpLocationScores, tmpAlingments, nTempAlignments);
                } else {
                    Log.Message("Skipping alignment computation.");
                }
                if (nTempAlignments > 0) {
                    verbose(0, "Aligned interval: ", tmpAlingments[nTempAlignments - 1].mappedInterval);
                }
            }
            if (pacbioDebug) {
                Log.Message("Alignment took %fs", tmr.ET());
            }

            read->AllocScores(tmpLocationScores, nTempAlignments);
            read->Alignments = tmpAlingments;

            delete[] tmpLocationScores;
            tmpLocationScores = 0;

            /**
             * Process alignments: removed short alignments, choos combination of aligned segments
             * that have the highest score and cover the read best
             */
            if (read->Calculated > 0) {
                bool mapped = reconcileRead(group);
                if (mapped) {
                    sortRead(group);
                } else {
                    verbose(0, true, "%s (%d) not mapped", read, read->length);
                }
                WriteRead(group->fullRead, mapped);
            } else {
                verbose(0, true, "%s (%d) not mapped", read, read->length);
                WriteRead(group->fullRead, false);
            }
        } else {
            verbose(0, true, "%s (%d) not mapped. No candidates found for read: unmapped.", read, read->length);
            //No candidates found for read
            WriteRead(group->fullRead, false);
        }

        delete[] anchorsFwd;
        anchorsFwd = 0;
        delete[] anchorsRev;
        anchorsRev = 0;

        if (intervals != 0) {
            for (int i = 0; i < nIntervals; ++i) {
                if (intervals[i] != 0) {
                    delete intervals[i];
                    intervals[i] = 0;
                }
            }
            delete[] intervals;
            intervals = 0;
            nIntervals = 0;
        }

        delete readCoordsTree;
        readCoordsTree = 0;
        processTime += overallTmr.ET();

        verbose(0, true, "###########################################");
        verbose(0, true, "###########################################");
        verbose(0, true, "###########################################");
        verbose(0, true, "");

    }
    this->longReadBatchIndex = 0;
    fprintf(stderr, "\tbatch processed YAY\n");
}

