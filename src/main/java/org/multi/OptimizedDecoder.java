package org.multi;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class OptimizedDecoder {

    
    
    
    private static final double[] POWERS_OF_10 = {
            1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0,
            10000000.0, 100000000.0, 1000000000.0, 10000000000.0,
            100000000000.0, 1000000000000.0, 10000000000000.0,
            100000000000000.0, 1000000000000000.0, 10000000000000000.0,
            100000000000000000.0, 1000000000000000000.0
    };

    public static Object[] decodeOptimized(String filePath) throws IOException {
        int totalDataPoints = 0;
        int totalMedoids = 0;

        
        
        
        int dim, packSize, blockSize, pageCount;
        List<Integer> pageKList = new ArrayList<>();
        List<Integer> pageDataPointsList = new ArrayList<>();
        List<int[]> maxDecimalList = new ArrayList<>();
        List<long[]> minValueList = new ArrayList<>();
        List<long[]> minBaseList = new ArrayList<>();
        List<int[]> maxBaseBitLenList = new ArrayList<>();
        List<int[][]> packResMetadataList = new ArrayList<>();

        BitStreamReader metaReader = new BitStreamReader(filePath);
        dim = (int) metaReader.readBits(8);
        packSize = (int) metaReader.readBits(16);
        blockSize = (int) metaReader.readBits(16);
        pageCount = (int) metaReader.readBits(32);

        for (int i = 0; i < pageCount; i++) {
            int k = (int) metaReader.readBits(16);
            int pageDataPoints = (int) metaReader.readBits(16);
            pageKList.add(k);
            pageDataPointsList.add(pageDataPoints);
            totalMedoids += k;
            totalDataPoints += pageDataPoints;

            int[] maxDecimals = new int[dim];
            long[] minValues = new long[dim];
            for (int j = 0; j < dim; j++) {
                maxDecimals[j] = (int) metaReader.readBits(8);
                int valBitLen = (int) metaReader.readBits(8);
                long sign = metaReader.readBits(1);
                long absValue = metaReader.readBits(valBitLen);
                minValues[j] = (sign == 0) ? absValue : -absValue;
            }
            maxDecimalList.add(maxDecimals);
            minValueList.add(minValues);

            long[] minBases = new long[dim];
            int[] maxBaseBits = new int[dim];
            for (int j = 0; j < dim; j++) {
                int valBitLen = (int) metaReader.readBits(8);
                long sign = metaReader.readBits(1);
                long absValue = metaReader.readBits(valBitLen);
                minBases[j] = (sign == 0) ? absValue : -absValue;
                maxBaseBits[j] = (int) metaReader.readBits(8);
            }
            minBaseList.add(minBases);
            maxBaseBitLenList.add(maxBaseBits);

            int packNum = (int) metaReader.readBits(16);
            int[][] packMetadata = new int[packNum][dim];
            for (int packIdx = 0; packIdx < packNum; packIdx++) {
                for (int j = 0; j < dim; j++) {
                    packMetadata[packIdx][j] = (int) metaReader.readBits(8);
                }
            }
            packResMetadataList.add(packMetadata);
        }

    
    
    
    double[][] finalDecodedData = new double[totalDataPoints][dim];
    double[][] finalMedoids = new double[totalMedoids][dim];
    long[] finalFrequencies = new long[totalMedoids];
    List<double[][]> finalMedoidsAsDoubleList = new ArrayList<>(pageCount);
    List<long[]> decodedFrequencyList = new ArrayList<>(pageCount);

    
    
    
    BitStreamReader dataReader = new BitStreamReader(filePath);
        
        dataReader.skipBits(8 + 16 + 16 + 32);
        for (int i = 0; i < pageCount; i++) {
            dataReader.skipBits(16 + 16);
            for (int j = 0; j < dim; j++) {
                dataReader.skipBits(8);
                int valBitLen = (int) dataReader.readBits(8);
                dataReader.skipBits(1 + valBitLen);
            }
            for (int j = 0; j < dim; j++) {
                int valBitLen = (int) dataReader.readBits(8);
                dataReader.skipBits(1 + valBitLen);
                dataReader.skipBits(8);
            }
            int packNum = (int) dataReader.readBits(16);
            dataReader.skipBits((long) packNum * dim * 8);
        }

        
        int dataCursor = 0;
        int medoidCursor = 0;

        
        for (int i = 0; i < pageCount; i++) {
            int k = pageKList.get(i);
            int pageDataPoints = pageDataPointsList.get(i);
            long[] minValues = minValueList.get(i);
            int[] maxDecimals = maxDecimalList.get(i);

            
            long[] minBases = minBaseList.get(i);
            int[] maxBaseBits = maxBaseBitLenList.get(i);
            double[][] currentPageDoubleMedoids = new double[k][dim];
            long[][] currentPageLongMedoids = new long[k][dim];

            for (int medoidIdx = 0; medoidIdx < k; medoidIdx++) {
                for (int j = 0; j < dim; j++) {
                    long offsetValue = dataReader.readBits(maxBaseBits[j]);
                    long longMedoidValue = offsetValue + minBases[j];
                    currentPageLongMedoids[medoidIdx][j] = longMedoidValue;

                    long integerValue = longMedoidValue + minValues[j];
                    double doubleValue = (double) integerValue / Math.pow(10, maxDecimals[j]);

                    
                    
                    
                    currentPageDoubleMedoids[medoidIdx][j] = doubleValue;
                    finalMedoids[medoidCursor + medoidIdx][j] = doubleValue;
                }
            }
            finalMedoidsAsDoubleList.add(currentPageDoubleMedoids);

            
            long[] clusterCounts = new long[k];
            if (k > 0) {
                int blockNum = (k + blockSize - 1) / blockSize;
                int[] blockMeta = new int[blockNum];
                for (int blockIdx = 0; blockIdx < blockNum; blockIdx++) {
                    blockMeta[blockIdx] = (int) dataReader.readBits(8);
                }
                for (int blockIdx = 0; blockIdx < blockNum; blockIdx++) {
                    int start = blockIdx * blockSize;
                    int end = Math.min(start + blockSize, k);
                    int bitsToRead = blockMeta[blockIdx];
                    for (int j = start; j < end; j++) {
                        clusterCounts[j] = dataReader.readBits(bitsToRead);
                        finalFrequencies[medoidCursor + j] = clusterCounts[j];
                    }
                }
            }
            long[] pageCumulativeFreq = new long[k];
            if (k > 0) {
                pageCumulativeFreq[0] = clusterCounts[0];
                for (int j = 1; j < k; j++) {
                    pageCumulativeFreq[j] = pageCumulativeFreq[j-1] + clusterCounts[j];
                }
            }
            decodedFrequencyList.add(pageCumulativeFreq);

            
            int[][] packMetadata = packResMetadataList.get(i);
            long[][] residualSeries = new long[pageDataPoints][dim];
            if (pageDataPoints > 0) {
                int dataCounter = 0;
                for (int packIdx = 0; packIdx < packMetadata.length; packIdx++) {
                    int[] bitsPerDim = packMetadata[packIdx];
                    int pointsInPack = (packIdx == packMetadata.length - 1)
                            ? (pageDataPoints - (packIdx * packSize))
                            : packSize;
                    for (int p = 0; p < pointsInPack; p++) {
                        for (int j = 0; j < dim; j++) {
                            long zigzagValue = dataReader.readBits(bitsPerDim[j]);
                            residualSeries[dataCounter][j] = residual.zigzagDecode(zigzagValue);
                        }
                        dataCounter++;
                    }
                }
            }

            
            int currentMedoidIndex = 0;
            for (int dataIdx = 0; dataIdx < pageDataPoints; dataIdx++) {
                if (k > 0 && currentMedoidIndex < k - 1 && dataIdx >= pageCumulativeFreq[currentMedoidIndex]) {
                    currentMedoidIndex++;
                }
                long[] basePoint = (k > 0) ? currentPageLongMedoids[currentMedoidIndex] : new long[dim];
                for (int j = 0; j < dim; j++) {
                    long integerValue = basePoint[j] + residualSeries[dataIdx][j] + minValues[j];
                    finalDecodedData[dataCursor + dataIdx][j] = (double) integerValue / Math.pow(10, maxDecimals[j]);
                }
            }

            
            dataCursor += pageDataPoints;
            medoidCursor += k;
        }

        
        
        
        return new Object[]{finalDecodedData, finalMedoidsAsDoubleList, decodedFrequencyList, finalMedoids, finalFrequencies};
    }

}
