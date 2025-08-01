package org.multi;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Encoder {
    public static Object[] generateBitstreamEncodedBasePoints(long[][] medoids, int dim) {
        BitBuffer baseBitstream = new BitBuffer();
        long[] minBase = new long[dim];
        int[] maxBaseBitLen = new int[dim];
        for(int i=0; i<dim; i++){
            minBase[i] =  medoids[0][i];
            maxBaseBitLen[i] = 0;
        }


        
        for (long[] point : medoids) {
            for(int i=0; i<dim; i++){
                if(point[i]<minBase[i]){
                    minBase[i] = point[i];
                }
            }
        }

        int length = medoids.length;
        for(int i=0; i<length; i++){
            for(int j=0; j<dim; j++){
                medoids[i][j] -= minBase[j];
                int mBit = support.bitsRequiredNoSign(medoids[i][j]);
                maxBaseBitLen[j] = Math.max(maxBaseBitLen[j], mBit);
            }
        }

        
        for (long[] point : medoids) {
            for (int i=0; i<dim; i++){
                support.appendToBitstream(baseBitstream, point[i], maxBaseBitLen[i]);
            }
        }

        return new Object[]{minBase, maxBaseBitLen, baseBitstream};
    }

    public static BitBuffer generateBitstreamFrequency(long[] clusterSize, int block_size){
        BitBuffer freq_bit = new BitBuffer();
        BitBuffer freq_meta = new BitBuffer();


        long[] clusterDelta = new long[clusterSize.length];
        clusterDelta[0] = clusterSize[0];
        for(int i=1; i<clusterSize.length; i++){
            clusterDelta[i] = clusterSize[i] - clusterSize[i-1];
        }




        int totalLength = clusterSize.length;
        int numBlocks = (totalLength + block_size - 1) / block_size;

        for (int blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
            int start = blockIdx * block_size;
            int end = Math.min(start + block_size, totalLength);

            long maxFrequency = 0;
            for (int i = start; i < end; i++) {
                long freq = clusterDelta[i];
                if (freq > maxFrequency) {
                    maxFrequency = freq;
                }
            }

            int maxFreqBit = support.bitsRequiredNoSign(maxFrequency);

            support.appendToBitstream(freq_meta, maxFreqBit, 8);

            for (int i = start; i < end; i++) {
                long freq = clusterDelta[i];
                support.appendToBitstream(freq_bit, freq, maxFreqBit);
            }
        }

        BitBuffer freq_res = new BitBuffer();
        freq_res.merge(freq_meta);
        freq_res.merge(freq_bit);
        return freq_res;
    }

    public static Object[] generateBitstreamEncodedDataWithHuffmanVaried(
            long[][] residualSeries, int[] clusterAssignment, String[] huffmanCode, int packSize, int dim) {
        BitBuffer residual_bit = new BitBuffer();
        BitBuffer base_bit = new BitBuffer();
        long totalResidualCost = 0;
        long totalBaseIdResCost = 0;
        int[][] packResMetadata = new int[(int) Math.ceil((double) residualSeries.length / packSize)][dim];
        int packIndex = 0;
        int totalDataPoints = residualSeries.length;
        int packNum = (totalDataPoints + packSize - 1) / packSize;

        for (int packStart = 0; packStart < residualSeries.length; packStart += packSize) {
            int packEnd = Math.min(packStart + packSize, residualSeries.length);
            int[] maxResidualBits = new int[dim];
            for(int i=0; i<dim; i++){
                maxResidualBits[i] = 0;
            }

            for (int i = packStart; i < packEnd; i++) {
                for(int j=0; j<dim; j++){
                    maxResidualBits[j] = Math.max(maxResidualBits[j], support.bitsRequiredNoSign(residualSeries[i][j]));
                }
            }

            for(int i=0; i<dim; i++){
                packResMetadata[packIndex][i] = maxResidualBits[i];
            }
            packIndex++;

            
            for (int i = packStart; i < packEnd; i++) {
                int baseId = clusterAssignment[i];
                String huffCode = huffmanCode[baseId];
                support.appendToBitstreamStr(base_bit, huffCode, huffCode.length());
                for(int j=0; j<dim; j++){
                    support.appendToBitstream(residual_bit, residualSeries[i][j], maxResidualBits[j]);
                    totalResidualCost += maxResidualBits[j];
                }
                totalBaseIdResCost += huffCode.length();
            }
        }

        return new Object[]{base_bit, residual_bit, packResMetadata, packNum, totalDataPoints};
    }

    public static Object[] saveCompressedDataWithBitstreamVaried(
            String filePath,
            long[][] medoids,
            int[] maxDecimalPlaces,
            int packSize,
            int blockSize,
            long[] minBase,
            int[] maxBaseBitLen,
            BitBuffer baseBitstream,
            BitBuffer HuffBitstream,
            BitBuffer baseresbit,
            BitBuffer residualbit,
            int[][] packResMetadata,
            int dim,
            int packNum,
            int totalDataPoints
    ) throws IOException {
        BitBuffer headerBitstream = new BitBuffer();
        BitBuffer packMetaBitstream = new BitBuffer();

        
        System.out.println("dim: " + dim);
        support.appendToBitstream(headerBitstream, dim, 8);
        System.out.println(headerBitstream.toBinaryString());
        for(int i=0; i<dim; i++){
            support.appendToBitstream(headerBitstream, maxDecimalPlaces[i], 8);
        }


        support.appendToBitstream(headerBitstream, medoids.length, 16);

        for(int i=0; i<dim; i++){
            int bMinBitLen = support.bitsRequiredNoSign(Math.abs(minBase[i])) + 1; 
            support.appendToBitstream(headerBitstream, bMinBitLen, 8);
            support.appendToBitstream(headerBitstream, minBase[i] >= 0 ? 0 : 1, 1);
            support.appendToBitstream(headerBitstream, Math.abs(minBase[i]), bMinBitLen);
        }

        for(int i=0; i<dim; i++){
            support.appendToBitstream(headerBitstream, maxBaseBitLen[i], 8);
        }

        support.appendToBitstream(packMetaBitstream, packSize, 16);
        support.appendToBitstream(packMetaBitstream, packNum, 32);
        support.appendToBitstream(packMetaBitstream, totalDataPoints, 32);
        for (int[] packMetadata : packResMetadata) {
            for(int i=0; i<dim; i++){
                support.appendToBitstream(packMetaBitstream, packMetadata[i], 8);
            }
        }

        
        BitBuffer finalBitstream = new BitBuffer();

        
        finalBitstream.merge(headerBitstream);
        finalBitstream.merge(baseBitstream);
        finalBitstream.merge(HuffBitstream);
        finalBitstream.merge(packMetaBitstream);
        finalBitstream.merge(baseresbit);
        finalBitstream.merge(residualbit);

        
        byte[] finalBytestream = finalBitstream.toByteArray();
        try (FileOutputStream fos = new FileOutputStream(filePath)) {
            fos.write(finalBytestream);
        }
        System.out.println("first byte: " + finalBytestream[0]);
        System.out.println(finalBytestream[1]);
        System.out.println(finalBytestream[2]);
        System.out.println(finalBytestream[3]);

        return new Object[]{finalBytestream.length, headerBitstream.size(), baseBitstream.size(), HuffBitstream.size(), packMetaBitstream.size(), baseresbit.size(), residualbit.size()};
    }

    public static Object[] clusterEncoder(String dataset, int k, int packSize, int blockSize, int dim) throws IOException {
        String inputFilePath = dataset + ".csv";
        String outputFilePath = dataset + "_compressed_varied.csv";
        String deltaOutputFile = dataset + "_delta.csv";
        String medoidsOutputFile = dataset + "_medoid.csv";
        String residualOutputFile = dataset + "_res.csv";
        Object[] preprocessingRes = Preprocessor.preprocessing(inputFilePath,dim);
        long[][] data = (long[][]) preprocessingRes[0];
        int[] maxDecimalPlacements = (int[]) preprocessingRes[1];
        System.out.println(Arrays.toString(maxDecimalPlacements));

        Object[] clusteringRes = Base.kMedoidLogCost(data,k,2,1.0e-3, dim);
        long[][] medoids = (long[][]) clusteringRes[0];
        int[] clusterAssignment = (int[]) clusteringRes[1];
        long[] clusterSize = (long[]) clusteringRes[2];




        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(medoidsOutputFile))) {
            for (long[] Point : medoids) {
                for(long aspect : Point){
                    writer.write(aspect + ", ");
                }
                writer.write("\n");  
            }
            writer.flush();  
        } catch (IOException e) {
            System.err.println("Error writing to file: " + e.getMessage());
        }

        System.out.println(Arrays.toString(clusterSize));
        BitBuffer frebit =generateBitstreamFrequency(clusterSize,blockSize);
        String[] canonicalCodes = Huff.generateHuffmanCodes(medoids,clusterSize);

        long[][] residualSeries = residual.residualCalculationZigzag(data,medoids,clusterAssignment,residualOutputFile,true,dim);







        Object[] bitstreamResult = generateBitstreamEncodedDataWithHuffmanVaried(residualSeries,clusterAssignment,canonicalCodes,packSize,dim);
        BitBuffer resbasebit = (BitBuffer) bitstreamResult[0];
        BitBuffer residualbit = (BitBuffer) bitstreamResult[1];
        int[][] packResMetadata = (int[][]) bitstreamResult[2];
        int packNum = (int) bitstreamResult[3];
        int totalDataPoints = (int) bitstreamResult[4];

        Object[] basePointResult = generateBitstreamEncodedBasePoints(medoids,dim);
        long[] minBase = (long[]) basePointResult[0];
        int[] maxBaseBitLen = (int[]) basePointResult[1];
        BitBuffer baseBitstream = (BitBuffer) basePointResult[2];

        return saveCompressedDataWithBitstreamVaried(outputFilePath,medoids,maxDecimalPlacements,packSize,blockSize,minBase,maxBaseBitLen,
                baseBitstream,frebit,resbasebit,residualbit,packResMetadata,dim, packNum, totalDataPoints);
    }

    public static Object[] clusterEncoderPage(String dataset, int k, int packSize, int blockSize, int dim, int pageSize, long[][] data, int[] maxDecimalPlacements) throws IOException {
        String inputFilePath = dataset + ".csv";
        String outputFilePath = dataset + "_compressed_paged.csv";
        String deltaOutputFile = dataset + "_delta.csv";
        String medoidsOutputFile = dataset + "_medoid.csv";
        String residualOutputFile = dataset + "_res.csv";

        





        
        int totalDataPoints = data.length;
        int totalPages = (totalDataPoints + pageSize - 1) / pageSize;
        System.out.println("Total pages: " + totalPages);

        
        List<long[][]> medoidsList = new ArrayList<>();
        List<BitBuffer> frebitList = new ArrayList<>();
        List<BitBuffer> resbasebitList = new ArrayList<>();
        List<BitBuffer> residualbitList = new ArrayList<>();
        List<int[][]> packResMetadataList = new ArrayList<>();
        List<long[]> minBaseList = new ArrayList<>();
        List<int[]> maxBaseBitLenList = new ArrayList<>();
        List<BitBuffer> baseBitstreamList = new ArrayList<>();
        List<Integer> pageDataPointsList = new ArrayList<>();
        List<Integer> pagePackNumList = new ArrayList<>();

        long limitT = 0;
        for(int i=0; i<dim; i++){
            limitT += 4;
        }
        
        for (int page = 0; page < totalPages; page++) {
            int start = page * pageSize;
            int end = Math.min(start + pageSize, totalDataPoints);
            
            long[][] pageData = new long[end - start][dim];
            for (int i = start; i < end; i++) {
                pageData[i - start] = data[i];
            }
            long partTime = System.nanoTime();
            
            Object[] clusteringRes = Base.kMedoidLogCost(pageData, k, 2, 1.0e-3, dim);






            long clusterTime = System.nanoTime();
            System.out.println("clustering time: " + (clusterTime - partTime));



            long[][] medoids = (long[][]) clusteringRes[0];
            int[] clusterAssignment = (int[]) clusteringRes[1];
            long[] clusterSize = (long[]) clusteringRes[2];


            










            partTime = System.nanoTime();
            
            BitBuffer frebit = generateBitstreamFrequency(clusterSize, blockSize);
            System.out.println("generateBitstream: " + (System.nanoTime()-partTime));

            
            partTime = System.nanoTime();
            String[] canonicalCodes = Huff.generateHuffmanCodes(medoids, clusterSize);
            System.out.println("huffgenerate: " + (System.nanoTime()-partTime));

            partTime = System.nanoTime();
            
            long[][] residualSeries = residual.residualCalculationZigzag(pageData, medoids, clusterAssignment, residualOutputFile, true, dim);
            System.out.println("residualzigzag: " + (System.nanoTime()-partTime));


            partTime = System.nanoTime();
            
            Object[] bitstreamResult = generateBitstreamEncodedDataWithHuffmanVaried(residualSeries, clusterAssignment, canonicalCodes, packSize, dim);
            System.out.println("encode: " + (System.nanoTime()-partTime));
            BitBuffer resbasebit = (BitBuffer) bitstreamResult[0];
            BitBuffer residualbit = (BitBuffer) bitstreamResult[1];
            int[][] packResMetadata = (int[][]) bitstreamResult[2];
            int pagePackNum = (int) bitstreamResult[3];
            int pageDataPoints = (int) bitstreamResult[4];

            
            Object[] basePointResult = generateBitstreamEncodedBasePoints(medoids, dim);
            long[] minBase = (long[]) basePointResult[0];
            int[] maxBaseBitLen = (int[]) basePointResult[1];
            BitBuffer baseBitstream = (BitBuffer) basePointResult[2];

            
            pageDataPointsList.add(pageDataPoints);
            pagePackNumList.add(pagePackNum);
            medoidsList.add(medoids);
            frebitList.add(frebit);
            resbasebitList.add(resbasebit);
            residualbitList.add(residualbit);
            packResMetadataList.add(packResMetadata);
            minBaseList.add(minBase);
            maxBaseBitLenList.add(maxBaseBitLen);
            baseBitstreamList.add(baseBitstream);

            System.out.println("Page " + (page + 1) + " processed. Data points: " + (end - start));
        }

        
        return saveCompressedDataWithBitstreamVariedPage(outputFilePath,
                medoidsList, maxDecimalPlacements, packSize, blockSize,
                minBaseList, maxBaseBitLenList, baseBitstreamList,
                frebitList, resbasebitList, residualbitList, packResMetadataList, dim,
                pageSize, totalPages, totalDataPoints, pageDataPointsList,
                pagePackNumList);
    }

    public static Object[] saveCompressedDataWithBitstreamVariedPage(
            String filePath,
            List<long[][]> medoidsList,
            int[] maxDecimalPlaces,
            int packSize,
            int blockSize,
            List<long[]> minBaseList,
            List<int[]> maxBaseBitLenList,
            List<BitBuffer> baseBitstreamList,
            List<BitBuffer> frebitList,
            List<BitBuffer> resbasebitList,
            List<BitBuffer> residualbitList,
            List<int[][]> packResMetadataList,
            int dim,
            int pageSize,
            int pageNum,
            int totalDataPoints,
            List<Integer> pageDataPointsList,
            List<Integer> pagePackNumList
    ) throws IOException{
        BitBuffer headerBitstream = new BitBuffer();
        BitBuffer packMetaBitstream = new BitBuffer();
        BitBuffer finalbaseBitstream = new BitBuffer();
        BitBuffer finalfrebit = new BitBuffer();
        BitBuffer finalresbasebit = new BitBuffer();
        BitBuffer finalresidualbit = new BitBuffer();

        
        
        System.out.println("page size: " + pageSize);
        System.out.println("page num: " + pageNum);
        System.out.println("dim: " + dim);
        support.appendToBitstream(headerBitstream, pageSize, 32);
        support.appendToBitstream(headerBitstream, pageNum, 32);
        support.appendToBitstream(headerBitstream, dim, 8);

        for(int i=0; i<dim; i++){
            support.appendToBitstream(headerBitstream, maxDecimalPlaces[i], 8);
        }


        support.appendToBitstream(packMetaBitstream, packSize, 16);
        support.appendToBitstream(packMetaBitstream, totalDataPoints, 32);
        for(int pageIndex=0; pageIndex<medoidsList.size(); pageIndex++){
            support.appendToBitstream(headerBitstream, medoidsList.get(pageIndex).length, 16);
            long[] minBase = minBaseList.get(pageIndex);
            for(int i=0; i<dim; i++){
                int bMinBitLen = support.bitsRequiredNoSign(Math.abs(minBase[i])) + 1; 
                support.appendToBitstream(headerBitstream, bMinBitLen, 8);
                support.appendToBitstream(headerBitstream, minBase[i] >= 0 ? 0 : 1, 1);
                support.appendToBitstream(headerBitstream, Math.abs(minBase[i]), bMinBitLen);
            }
            int[] maxBaseBitLen = maxBaseBitLenList.get(pageIndex);
            for(int i=0; i<dim; i++){
                support.appendToBitstream(headerBitstream, maxBaseBitLen[i], 8);
            }
            int[][] packResMetadata = packResMetadataList.get(pageIndex);
            for (int[] packMetadata : packResMetadata) {
                for(int i=0; i<dim; i++){
                    support.appendToBitstream(packMetaBitstream, packMetadata[i], 8);
                }
            }

            finalbaseBitstream.merge(baseBitstreamList.get(pageIndex));
            finalfrebit.merge(frebitList.get(pageIndex));
            finalresbasebit.merge(resbasebitList.get(pageIndex));
            finalresidualbit.merge(residualbitList.get(pageIndex));
        }
        
        BitBuffer finalBitstream = new BitBuffer();

        
        finalBitstream.merge(headerBitstream);
        finalBitstream.merge(finalbaseBitstream);
        finalBitstream.merge(finalfrebit);
        finalBitstream.merge(packMetaBitstream);
        finalBitstream.merge(finalresbasebit);
        finalBitstream.merge(finalresidualbit);

        
        byte[] finalBytestream = finalBitstream.toByteArray();
        try (FileOutputStream fos = new FileOutputStream(filePath)) {
            fos.write(finalBytestream);
        }
        System.out.println("first byte: " + finalBytestream[0]);
        System.out.println(finalBytestream[1]);
        System.out.println(finalBytestream[2]);
        System.out.println(finalBytestream[3]);

        return new Object[]{finalBytestream.length, headerBitstream.size(), finalbaseBitstream.size(),
                finalfrebit.size(), packMetaBitstream.size(), finalresbasebit.size(), finalresidualbit.size()};
    }

    private static void saveMedoidsToCSV(long[][] medoids, String filePath) throws IOException {
        
        try (PrintWriter writer = new PrintWriter(new FileWriter(filePath))) {
            
            for (long[] row : medoids) {
                
                
                
                
                String line = Arrays.stream(row)
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(","));

                
                writer.println(line);
            }
        }
        
        System.out.println("✅ Medoids data successfully saved to: " + filePath);
    }


    public static Object[] clusterEncoderPageNoHuff(double[][] raw_data, int k, int packSize, int dim, int blockSize, String filePath, int method_idx) throws IOException {


        
        Object[] preprocessingRes = Preprocessor.preprocessingPos(raw_data, dim);
        long[][] data = (long[][]) preprocessingRes[0];
        int[] maxDecimalPlacements = (int[]) preprocessingRes[1];
        long[] minValues = (long[]) preprocessingRes[2];
        int numRows = (int) preprocessingRes[3];







        
        Object[] clusteringRes = new Object[3];
        if(method_idx==0){
            clusteringRes = Base.kMedoidLogCost_SIMD(data,k,2,1.0e-3,dim);
        }else if(method_idx==1){
            clusteringRes = Base.adaptiveGreedyBasePointSelection_simd(data,dim);
        }else if(method_idx==2){
            clusteringRes = oneDSolver.optimalKMedoid1D(data,k,1);
        }else if(method_idx==3){
            clusteringRes = oneDSolver.optimalAdaptiveBPRE1D(data,1);
        }






















        long[][] medoids = (long[][]) clusteringRes[0];
        int[] clusterAssignment = (int[]) clusteringRes[1];


        long[] clusterSize = (long[]) clusteringRes[2];













        int page_k = medoids.length;
        
        BitBuffer frebit = generateBitstreamFrequency(clusterSize, blockSize);

        
        String residualOutputFile = filePath + "_res.csv";

        long[][] residualSeries = residual.residualCalculationZigzagNoHuff_sorted(data, medoids, clusterAssignment, clusterSize, residualOutputFile, false, dim);

        Object[] bitstreamResult = generateBitstreamEncodedData(residualSeries, packSize, dim);


        BitBuffer resbit = (BitBuffer) bitstreamResult[0];
        int[][] packResMetadata = (int[][]) bitstreamResult[1];
        int pagePackNum = (int) bitstreamResult[2];
        int pageDataPoints = (int) bitstreamResult[3];

        
        Object[] basePointResult = generateBitstreamEncodedBasePointsOptimized(medoids, dim);
        long[] minBase = (long[]) basePointResult[0];
        int[] maxBaseBitLen = (int[]) basePointResult[1];
        BitBuffer baseBitstream = (BitBuffer) basePointResult[2];

        
        return new Object[]{maxDecimalPlacements, minValues, frebit, resbit, packResMetadata, pagePackNum, pageDataPoints, minBase, maxBaseBitLen, baseBitstream, page_k};
    }

    public static Object[] generateBitstreamEncodedBasePointsOptimized(long[][] medoids, int dim) {
        if (medoids == null || medoids.length == 0) {
            return new Object[]{new long[dim], new int[dim], new BitBuffer()};
        }

        int k = medoids.length;
        long[] minBase = new long[dim];
        Arrays.fill(minBase, Long.MAX_VALUE);

        for (long[] point : medoids) {
            for (int i = 0; i < dim; i++) {
                minBase[i] = Math.min(minBase[i], point[i]);
            }
        }

        long[][] offsetMedoids = new long[k][dim];
        int[] maxBaseBitLen = new int[dim];

        for (int i = 0; i < k; i++) {
            for (int j = 0; j < dim; j++) {
                long offsetValue = medoids[i][j] - minBase[j];
                offsetMedoids[i][j] = offsetValue;

                maxBaseBitLen[j] = Math.max(maxBaseBitLen[j], support.bitsRequiredNoSign(offsetValue));
            }
        }

        BitBuffer baseBitstream = new BitBuffer();
        for (long[] offsetPoint : offsetMedoids) {
            for (int i = 0; i < dim; i++) {
                support.appendToBitstream(baseBitstream, offsetPoint[i], maxBaseBitLen[i]);
            }
        }

        return new Object[]{minBase, maxBaseBitLen, baseBitstream};
    }

    public static Object[] generateBitstreamEncodedData(long[][] residualSeries, int packSize, int dim) {
        int totalDataPoints = residualSeries.length;
        if (totalDataPoints == 0) {
            return new Object[]{new BitBuffer(), new int[0][dim], 0, 0};
        }

        
        int packNum = (totalDataPoints + packSize - 1) / packSize;

        
        BitBuffer residual_bit = new BitBuffer();
        int[][] packResMetadata = new int[packNum][dim];

        int packIndex = 0;
        
        for (int packStart = 0; packStart < totalDataPoints; packStart += packSize) {
            int packEnd = Math.min(packStart + packSize, totalDataPoints);

            
            int[] maxResidualBits = new int[dim]; 
            for (int i = packStart; i < packEnd; i++) {
                for (int j = 0; j < dim; j++) {
                    int bitsNeeded = support.bitsRequiredNoSign(residualSeries[i][j]);
                    maxResidualBits[j] = Math.max(maxResidualBits[j], bitsNeeded);
                }
            }

            
            packResMetadata[packIndex] = maxResidualBits;

            for (int i = packStart; i < packEnd; i++) {
                for (int j = 0; j < dim; j++) {
                    support.appendToBitstream(residual_bit,residualSeries[i][j], maxResidualBits[j]);
                }
            }

            packIndex++;
        }

        return new Object[]{residual_bit, packResMetadata, packNum, totalDataPoints};
    }

    public static Object[] saveCompressedDataNoHuff(
            String outputFilePath, int dim, int packSize, int blockSize, int pageCount,
            List<int[]> maxDecimalList, List<long[]> minValueList, List<Integer> pageDataPointsList,
            List<BitBuffer> frebitList, List<BitBuffer> residualbitList, List<int[][]> packResMetadataList,
            List<long[]> minBaseList, List<int[]> maxBaseBitLenList, List<BitBuffer> baseBitstreamList,
            List<Integer> pageKList) throws IOException {

        
        BitBuffer headerBitstream = new BitBuffer();

        BitBuffer finalbaseBitstream = new BitBuffer();
        BitBuffer finalfrebit = new BitBuffer();
        BitBuffer finalresidualbit = new BitBuffer();

        support.appendToBitstream(headerBitstream, dim, 8);          
        support.appendToBitstream(headerBitstream, packSize, 16);    
        support.appendToBitstream(headerBitstream, blockSize, 16);   
        support.appendToBitstream(headerBitstream, pageCount-1, 32);   
        
        for (int i = 0; i < pageCount-1; i++) {
            
            int k = pageKList.get(i);
            support.appendToBitstream(headerBitstream, (long)k, 16); 
            support.appendToBitstream(headerBitstream, pageDataPointsList.get(i), 16); 

            
            int[] maxDecimals = maxDecimalList.get(i);
            long[] minValues = minValueList.get(i);
            for (int j = 0; j < dim; j++) {
                support.appendToBitstream(headerBitstream, maxDecimals[j], 8); 

                int minValueBitLen = support.bitsRequiredNoSign(Math.abs(minValues[j])) + 1;

                support.appendToBitstream(headerBitstream, minValueBitLen, 8);
                support.appendToBitstream(headerBitstream, minValues[j] >= 0 ? 0 : 1, 1);
                support.appendToBitstream(headerBitstream, Math.abs(minValues[j]), minValueBitLen);
            }

            
            long[] minBases = minBaseList.get(i);

            int[] maxBaseBits = maxBaseBitLenList.get(i);
            for (int j = 0; j < dim; j++) {
                int minBaseBitLen = support.bitsRequiredNoSign(Math.abs(minBases[j])) + 1;
                support.appendToBitstream(headerBitstream, minBaseBitLen, 8);
                support.appendToBitstream(headerBitstream, minBases[j] >= 0 ? 0 : 1, 1);
                support.appendToBitstream(headerBitstream, Math.abs(minBases[j]), minBaseBitLen);
                support.appendToBitstream(headerBitstream, maxBaseBits[j], 8);     
            }

            
            int[][] packMetadata = packResMetadataList.get(i);
            
            support.appendToBitstream(headerBitstream, packMetadata.length, 16);
            for (int[] pack : packMetadata) {
                for (int j = 0; j < dim; j++) {
                    support.appendToBitstream(headerBitstream, pack[j], 8); 
                }
            }

            finalfrebit.merge(frebitList.get(i));

            finalresidualbit.merge(residualbitList.get(i));

            finalbaseBitstream.merge(baseBitstreamList.get(i));
        }


        
        BitBuffer finalOutputBitstream = new BitBuffer();
        finalOutputBitstream.merge(headerBitstream);        
        finalOutputBitstream.merge(finalbaseBitstream);     
        finalOutputBitstream.merge(finalfrebit);   
        finalOutputBitstream.merge(finalresidualbit); 

        
        byte[] finalBytestream = finalOutputBitstream.toByteArray();
        try (FileOutputStream fos = new FileOutputStream(outputFilePath)) {
            fos.write(finalBytestream);
        }

        
        return new Object[]{
                (long) finalBytestream.length,      
                headerBitstream.size(),             
                finalbaseBitstream.size(),          
                finalfrebit.size(),        
                finalresidualbit.size()       
        };
    }

    public static Object[] saveCompressedDataNoHuff_Page(
            String outputFilePath, int dim, int packSize, int blockSize, int pageCount,
            List<int[]> maxDecimalList, List<long[]> minValueList, List<Integer> pageDataPointsList,
            List<BitBuffer> frebitList, List<BitBuffer> residualbitList, List<int[][]> packResMetadataList,
            List<long[]> minBaseList, List<int[]> maxBaseBitLenList, List<BitBuffer> baseBitstreamList,
            List<Integer> pageKList, int PageSize) throws IOException {

        
        BitBuffer finalOutputBitstream = new BitBuffer();

        
        long totalHeaderSize = 0;
        long totalBaseSize = 0;
        long totalFreSize = 0;
        long totalResidualSize = 0;

        
        BitBuffer globalHeader = new BitBuffer();
        support.appendToBitstream(globalHeader, dim, 8);          
        support.appendToBitstream(globalHeader, packSize, 16);    
        support.appendToBitstream(globalHeader, blockSize, 16);   
        support.appendToBitstream(globalHeader, pageCount, 32);   
        support.appendToBitstream(globalHeader, PageSize, 16);

        
        finalOutputBitstream.merge(globalHeader);
        totalHeaderSize += globalHeader.size();


        
        for (int i = 0; i < pageCount; i++) {
            
            BitBuffer pageMetadataBitstream = new BitBuffer();

            int k = pageKList.get(i);

            support.appendToBitstream(pageMetadataBitstream, (long)k, 16);
            support.appendToBitstream(pageMetadataBitstream, pageDataPointsList.get(i), 16);

            int[] maxDecimals = maxDecimalList.get(i);
            long[] minValues = minValueList.get(i);
            for (int j = 0; j < dim; j++) {
                support.appendToBitstream(pageMetadataBitstream, maxDecimals[j], 8);
                int minValueBitLen = support.bitsRequiredNoSign(Math.abs(minValues[j])) + 1;
                support.appendToBitstream(pageMetadataBitstream, minValueBitLen, 8);
                support.appendToBitstream(pageMetadataBitstream, minValues[j] >= 0 ? 0 : 1, 1);
                support.appendToBitstream(pageMetadataBitstream, Math.abs(minValues[j]), minValueBitLen);
            }

            long[] minBases = minBaseList.get(i);
            int[] maxBaseBits = maxBaseBitLenList.get(i);
            for (int j = 0; j < dim; j++) {
                int minBaseBitLen = support.bitsRequiredNoSign(Math.abs(minBases[j])) + 1;
                support.appendToBitstream(pageMetadataBitstream, minBaseBitLen, 8);
                support.appendToBitstream(pageMetadataBitstream, minBases[j] >= 0 ? 0 : 1, 1);
                support.appendToBitstream(pageMetadataBitstream, Math.abs(minBases[j]), minBaseBitLen);
                support.appendToBitstream(pageMetadataBitstream, maxBaseBits[j], 8);
            }

            int[][] packMetadata = packResMetadataList.get(i);
            support.appendToBitstream(pageMetadataBitstream, packMetadata.length, 16);
            for (int[] pack : packMetadata) {
                for (int j = 0; j < dim; j++) {
                    support.appendToBitstream(pageMetadataBitstream, pack[j], 8);
                }
            }

            
            BitBuffer baseData = baseBitstreamList.get(i);
            BitBuffer freData = frebitList.get(i);
            BitBuffer resData = residualbitList.get(i);

            
            finalOutputBitstream.merge(pageMetadataBitstream);
            finalOutputBitstream.merge(baseData);
            finalOutputBitstream.merge(freData);
            finalOutputBitstream.merge(resData);

            
            totalHeaderSize += pageMetadataBitstream.size();
            totalBaseSize += baseData.size();
            totalFreSize += freData.size();
            totalResidualSize += resData.size();
        }

        
        byte[] finalBytestream = finalOutputBitstream.toByteArray();
        try (FileOutputStream fos = new FileOutputStream(outputFilePath)) {
            fos.write(finalBytestream);
        }

        
        return new Object[]{
                (long) finalBytestream.length,
                totalHeaderSize,
                totalBaseSize,
                totalFreSize,
                totalResidualSize
        };
    }
}
