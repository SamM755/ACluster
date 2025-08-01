package org.multi;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.ByteArrayOutputStream;
import java.util.List;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorSpecies;
import jdk.incubator.vector.VectorOperators;

public class ClusterDecoder {
    public static double[][] convertLongToDouble(long[][] longValues, long[] maxDecimalPlaces) {
        if (longValues.length == 0 || longValues[0].length != maxDecimalPlaces.length) {
            throw new IllegalArgumentException("Input dimensions do not match maxDecimalPlaces length.");
        }

        int rows = longValues.length;
        int dim = maxDecimalPlaces.length;
        double[][] doubleValues = new double[rows][dim];

        double[] scaleFactors = new double[dim];
        for (int j = 0; j < dim; j++) {
            scaleFactors[j] = Math.pow(10, maxDecimalPlaces[j]); 
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < dim; j++) {
                doubleValues[i][j] = longValues[i][j] / scaleFactors[j]; 
            }
        }
        return doubleValues;
    }

    public static Object[] decodeBases(String filePath) throws IOException {
        BitStreamReader reader = new BitStreamReader(filePath);

        long dim = reader.readBits(8);

        long[] maxDecimalPlaces = new long[(int) dim];
        for(int i=0; i<dim; i++){
            maxDecimalPlaces[i] = reader.readBits(8);
        }

        int k = (int) reader.readBits(16);
        long[] minBase = new long[(int) dim];
        int[] bMinBitLen = new int[(int) dim];
        for(int i=0; i<dim; i++){
            bMinBitLen[i] = (int) reader.readBits(8);
            long sign = reader.readBits(1);
            long tmp_v = reader.readBits(bMinBitLen[i]);
            if(sign == 0){
                minBase[i] = tmp_v;
            }else{
                minBase[i] = -1 * tmp_v;
            }
        }
        int[] maxBaseBitLen = new int[(int)dim];
        for(int i=0; i<dim; i++){
            maxBaseBitLen[i] = (int) reader.readBits(8);
        }

        long[][] medoids = new long[k][(int)dim];
        for(int i=0; i<k; i++){
            for(int j=0; j<dim; j++){
                long tmp = reader.readBits(maxBaseBitLen[j]);
                medoids[i][j] = tmp + minBase[j];
            }
        }

        long[] frequency = new long[k];
        int blockSize = (int) reader.readBits(8);
        int blockNum = (k + blockSize - 1) / blockSize;
        int[] blockMeta = new int[blockNum];
        for(int i=0; i<blockNum; i++){
            blockMeta[i] = (int)reader.readBits(8);
        }
        for(int i=0; i<blockNum; i++){
            int start = i * blockSize;
            int end = Math.min(start + blockSize, k);
            for(int j=start; j<end; j++){
                frequency[j] = reader.readBits(blockMeta[i]);
            }
        }


        for(int i=1; i<k; i++){
            frequency[i] += frequency[i-1];
        }

        return new Object[]{medoids, frequency, maxDecimalPlaces, dim};
    }

    public static Object[] decodeAll(String filePath) throws IOException {
        BitStreamReader reader = new BitStreamReader(filePath);

        long dim = reader.readBits(8);

        long[] maxDecimalPlaces = new long[(int) dim];
        for(int i=0; i<dim; i++){
            maxDecimalPlaces[i] = reader.readBits(8);

        }

        int k = (int) reader.readBits(16);
        long[] minBase = new long[(int) dim];
        int[] bMinBitLen = new int[(int) dim];
        for(int i=0; i<dim; i++){
            bMinBitLen[i] = (int) reader.readBits(8);
            long sign = reader.readBits(1);
            long tmp_v = reader.readBits(bMinBitLen[i]);
            if(sign == 0){
                minBase[i] = tmp_v;
            }else{
                minBase[i] = -1 * tmp_v;
            }
        }
        int[] maxBaseBitLen = new int[(int)dim];
        for(int i=0; i<dim; i++){
            maxBaseBitLen[i] = (int) reader.readBits(8);
        }


        long[][] medoids = new long[k][(int)dim];
        for(int i=0; i<k; i++){
            for(int j=0; j<dim; j++){
                long tmp = reader.readBits(maxBaseBitLen[j]);
                medoids[i][j] = tmp + minBase[j];
            }
        }


        long[] frequency = new long[k];
        int blockSize = (int) reader.readBits(8);
        int blockNum = (k + blockSize - 1) / blockSize;
        int[] blockMeta = new int[blockNum];
        for(int i=0; i<blockNum; i++){
            blockMeta[i] = (int)reader.readBits(8);
        }
        for(int i=0; i<blockNum; i++){
            int start = i * blockSize;
            int end = Math.min(start + blockSize, k);
            for(int j=start; j<end; j++){
                frequency[j] = reader.readBits(blockMeta[i]);
            }
        }


        for(int i=1; i<k; i++){
            frequency[i] += frequency[i-1];
        }

        String[] huffmanCodes = Huff.generateHuffmanCodes(medoids,frequency);

        
        int packSize = (int) reader.readBits(16);
        int packNum = (int) reader.readBits(32);
        int totalDataPoints = (int) reader.readBits(32);
        int[][] packMetadata = new int[packNum][(int)dim];
        for(int packIndex=0; packIndex<packNum; packIndex++){
            for(int i=0; i<dim; i++){
                packMetadata[packIndex][i] = (int) reader.readBits(8);
            }
        }

        
        int[] index = new int[totalDataPoints];
        int counter = 0;
        StringBuilder currentBitstream = new StringBuilder();
        while (counter<totalDataPoints) {
            long nextBit = reader.readBits(1);
            currentBitstream.append(nextBit == 0 ? "0" : "1");

            
            for (int i = 0; i < huffmanCodes.length; i++) {
                if (currentBitstream.length() >= huffmanCodes[i].length() &&
                        currentBitstream.toString().startsWith(huffmanCodes[i])) {
                    
                    index[counter] = i;  
                    
                    currentBitstream.delete(0, huffmanCodes[i].length());
                    counter++;
                    break;
                }
            }
        }

        
        int data_cnt = 0;
        long[][] zigzagedSeries = new long[totalDataPoints][(int)dim];
        for(int packIdx = 0; packIdx<packNum; packIdx++){
            int packEnd = Math.min(totalDataPoints - (packIdx * packSize), packSize);
            int[] maxResidualBits = new int[(int)dim];
            for(int i=0; i<dim; i++){
                maxResidualBits[i] = packMetadata[packIdx][i];
            }
            for(int startIdx = 0; startIdx<packEnd; startIdx++){
                for(int i=0; i<dim; i++){
                    zigzagedSeries[data_cnt][i] = reader.readBits(maxResidualBits[i]);
                }
                data_cnt++;
            }
        }
        long[][] residualSeries = new long[totalDataPoints][(int)dim];
        long[][] decodedData = new long[totalDataPoints][(int)dim];
        for(int j=0; j<totalDataPoints; j++){
            long[] base = medoids[index[j]];
            for(int i=0; i<dim; i++){
                residualSeries[j][i] = residual.zigzagDecode(zigzagedSeries[j][i]);
                decodedData[j][i] = base[i] + residualSeries[j][i];
            }

        }


        return new Object[]{decodedData, medoids, frequency, maxDecimalPlaces, dim, index};
    }

    public static Object[] decodeMedoidPage(String filePath) throws IOException{
        BitStreamReader reader = new BitStreamReader(filePath);
        long pageSize = reader.readBits(32);
        long pageNum = reader.readBits(32);
        long dim = reader.readBits(8);



        long[] maxDecimalPlaces = new long[(int) dim];
        for(int i=0; i<dim; i++){
            maxDecimalPlaces[i] = reader.readBits(8);

        }


        List<long[][]> medoidsList = new ArrayList<>();
        List<long[]> minBaseList= new ArrayList<>();
        List<int[]> maxBaseBitLenList= new ArrayList<>();
        List<long[]> frequencyList= new ArrayList<>();
        List<String[]> huffmanCodesList = new ArrayList<>();
        List<Integer> kPage = new ArrayList<>();
        List<int[][]> packResMetadataList= new ArrayList<>();



        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = (int) reader.readBits(16);

            kPage.add(k);
            long[] minBase = new long[(int) dim];
            int[] bMinBitLen = new int[(int) dim];
            for(int i=0; i<dim; i++){
                bMinBitLen[i] = (int) reader.readBits(8);
                long sign = reader.readBits(1);
                long tmp_v = reader.readBits(bMinBitLen[i]);
                if(sign == 0){
                    minBase[i] = tmp_v;
                }else{
                    minBase[i] = -1 * tmp_v;
                }
            }
            minBaseList.add(minBase);
            int[] maxBaseBitLen = new int[(int)dim];
            for(int i=0; i<dim; i++){
                maxBaseBitLen[i] = (int) reader.readBits(8);
            }


            maxBaseBitLenList.add(maxBaseBitLen);
        }
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = kPage.get(pageIndex);
            long[] minBase_page = minBaseList.get(pageIndex);
            int[] maxBaseBitLen_page = maxBaseBitLenList.get(pageIndex);
            long[][] medoids = new long[k][(int)dim];
            for(int i=0; i<k; i++){
                for(int j=0; j<dim; j++){
                    long tmp = reader.readBits(maxBaseBitLen_page[j]);
                    medoids[i][j] = tmp + minBase_page[j];
                }
            }


            medoidsList.add(medoids);
        }
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = kPage.get(pageIndex);
            long[] frequency = new long[k];
            int blockSize = (int) reader.readBits(8);

            int blockNum = (k + blockSize - 1) / blockSize;
            int[] blockMeta = new int[blockNum];
            for(int i=0; i<blockNum; i++){
                blockMeta[i] = (int)reader.readBits(8);
            }

            for(int i=0; i<blockNum; i++){
                int start = i * blockSize;
                int end = Math.min(start + blockSize, k);
                for(int j=start; j<end; j++){
                    frequency[j] = reader.readBits(blockMeta[i]);
                }
            }
            for(int i=1; i<k; i++){
                frequency[i] += frequency[i-1];
            }


            String[] huffmanCodes = Huff.generateHuffmanCodes(medoidsList.get(pageIndex), frequency);
            huffmanCodesList.add(huffmanCodes);
            frequencyList.add(frequency);
        }
        int k_cnt = 0;
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            k_cnt += kPage.get(pageIndex);
        }
        long[][] full_medoid_list = new long[k_cnt][(int)dim];
        long[] full_frequency = new long[k_cnt];
        int med_cnt = 0;
        for(int i=0; i<pageNum; i++){
            long[] cur_frequency = frequencyList.get(i);
            long[][] medoids = medoidsList.get(i);
            for(int j=0; j<cur_frequency.length; j++){
                full_medoid_list[med_cnt] = Arrays.copyOf(medoids[j], medoids[j].length);
                full_frequency[med_cnt] = cur_frequency[j];
                med_cnt++;
            }
        }
        return new Object[]{full_medoid_list, full_frequency, maxDecimalPlaces, dim};
    }

    public static Object[] decodeAllPage(String filePath) throws IOException {
        BitStreamReader reader = new BitStreamReader(filePath);
        long pageSize = reader.readBits(32);
        long pageNum = reader.readBits(32);
        long dim = reader.readBits(8);



        long[] maxDecimalPlaces = new long[(int) dim];
        for(int i=0; i<dim; i++){
            maxDecimalPlaces[i] = reader.readBits(8);

        }


        List<long[][]> medoidsList = new ArrayList<>();
        List<long[]> minBaseList= new ArrayList<>();
        List<int[]> maxBaseBitLenList= new ArrayList<>();
        List<long[]> frequencyList= new ArrayList<>();
        List<String[]> huffmanCodesList = new ArrayList<>();

        List<int[][]> packResMetadataList= new ArrayList<>();

        List<Integer> kList = new ArrayList<>();

        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = (int) reader.readBits(16);
            long[] minBase = new long[(int) dim];
            int[] bMinBitLen = new int[(int) dim];
            for(int i=0; i<dim; i++){
                bMinBitLen[i] = (int) reader.readBits(8);
                System.out.println(bMinBitLen[i]);
                long sign = reader.readBits(1);
                long tmp_v = reader.readBits(bMinBitLen[i]);
                if(sign == 0){
                    minBase[i] = tmp_v;
                }else{
                    minBase[i] = -1 * tmp_v;
                }
            }
            minBaseList.add(minBase);
            int[] maxBaseBitLen = new int[(int)dim];
            for(int i=0; i<dim; i++){
                maxBaseBitLen[i] = (int) reader.readBits(8);
            }


            maxBaseBitLenList.add(maxBaseBitLen);
            kList.add(k);
        }
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = kList.get(pageIndex);
            long[] minBase_page = minBaseList.get(pageIndex);
            int[] maxBaseBitLen_page = maxBaseBitLenList.get(pageIndex);
            long[][] medoids = new long[k][(int)dim];
            for(int i=0; i<k; i++){
                for(int j=0; j<dim; j++){
                    long tmp = reader.readBits(maxBaseBitLen_page[j]);
                    medoids[i][j] = tmp + minBase_page[j];
                }
            }


            medoidsList.add(medoids);
        }
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int k = kList.get(pageIndex);
            long[] frequency = new long[k];
            int blockSize = (int) reader.readBits(8);

            int blockNum = (k + blockSize - 1) / blockSize;
            int[] blockMeta = new int[blockNum];
            for(int i=0; i<blockNum; i++){
                blockMeta[i] = (int)reader.readBits(8);
            }

            for(int i=0; i<blockNum; i++){
                int start = i * blockSize;
                int end = Math.min(start + blockSize, k);
                for(int j=start; j<end; j++){
                    frequency[j] = reader.readBits(blockMeta[i]);
                }
            }
            for(int i=1; i<k; i++){
                frequency[i] += frequency[i-1];
            }


            String[] huffmanCodes = Huff.generateHuffmanCodes(medoidsList.get(pageIndex), frequency);
            huffmanCodesList.add(huffmanCodes);
            frequencyList.add(frequency);
        }
        int packSize = (int) reader.readBits(16);

        int totalDataPoints = (int) reader.readBits(32);

        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            long actualPageSize = Math.min(pageSize, totalDataPoints - pageIndex * pageSize);;

            long packNum = (actualPageSize + packSize - 1) / packSize;

            int[][] packResMeta_page = new int[(int)packNum][(int)dim];
            for(int packIndex=0; packIndex<packNum; packIndex++){
                for(int i=0; i<dim; i++){
                    packResMeta_page[packIndex][i] = (int) reader.readBits(8);
                }
            }
            packResMetadataList.add(packResMeta_page);

        }
        List<int[]> medoidIndexList= new ArrayList<>();
        int k_start = 0;
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            String[] huffmanCodes = huffmanCodesList.get(pageIndex);
            long actualPageSize = Math.min(pageSize, totalDataPoints - pageIndex * pageSize);;
            int[] index = new int[(int)actualPageSize];
            int counter = 0;
            StringBuilder currentBitstream = new StringBuilder();
            while (counter<actualPageSize) {
                long nextBit = reader.readBits(1);
                currentBitstream.append(nextBit == 0 ? "0" : "1");

                
                for (int i = 0; i < huffmanCodes.length; i++) {
                    if (currentBitstream.length() >= huffmanCodes[i].length() &&
                            currentBitstream.toString().startsWith(huffmanCodes[i])) {
                        
                        index[counter] = i + k_start;  
                        
                        currentBitstream.delete(0, huffmanCodes[i].length());
                        counter++;
                        break;
                    }
                }
            }
            medoidIndexList.add(index);
            k_start += kList.get(pageIndex);

        }
        int full_data_cnt = 0;
        long[][] full_decoded_data = new long[totalDataPoints][(int)dim];
        long[] full_index_list = new long[totalDataPoints];
        k_start = 0;
        for(int pageIndex=0; pageIndex<pageNum; pageIndex++){
            int[] index = medoidIndexList.get(pageIndex);
            int[][] packMetadata = packResMetadataList.get(pageIndex);
            long[][] medoids = medoidsList.get(pageIndex);
            long actualPageSize = Math.min(pageSize, totalDataPoints - pageIndex * pageSize);
            long packNum = (actualPageSize + packSize - 1) / packSize;

            int data_cnt = 0;
            long[][] zigzagedSeries = new long[(int)actualPageSize][(int)dim];
            for(int packIdx = 0; packIdx<packNum; packIdx++){
                int packEnd = Math.min((int)actualPageSize - (packIdx * packSize), packSize);
                int[] maxResidualBits = new int[(int)dim];
                for(int i=0; i<dim; i++){
                    maxResidualBits[i] = packMetadata[packIdx][i];
                }
                for(int startIdx = 0; startIdx<packEnd; startIdx++){
                    for(int i=0; i<dim; i++){
                        zigzagedSeries[data_cnt][i] = reader.readBits(maxResidualBits[i]);
                    }
                    data_cnt++;
                }
            }


            long[][] residualSeries = new long[(int)actualPageSize][(int)dim];
            long[][] decodedData = new long[(int)actualPageSize][(int)dim];
            for(int j=0; j<actualPageSize; j++){
                long[] base = medoids[index[j]-k_start];
                for(int i=0; i<dim; i++){
                    residualSeries[j][i] = residual.zigzagDecode(zigzagedSeries[j][i]);
                    decodedData[j][i] = base[i] + residualSeries[j][i];
                }
                full_decoded_data[full_data_cnt] = Arrays.copyOf(decodedData[j], decodedData[j].length);
                full_index_list[full_data_cnt] = index[j];
                full_data_cnt++;
            }
            k_start += kList.get(pageIndex);
        }
        int total_k = 0;
        for(int i=0; i<kList.size(); i++){
            total_k += kList.get(i);
        }
        long[][] full_medoid_list = new long[total_k][(int)dim];
        long[] full_frequency = new long[total_k];
        int med_cnt = 0;
        for(int i=0; i<pageNum; i++){
            long[] cur_frequency = frequencyList.get(i);
            long[][] medoids = medoidsList.get(i);
            for(int j=0; j<cur_frequency.length; j++){
                full_medoid_list[med_cnt] = Arrays.copyOf(medoids[j], medoids[j].length);
                full_frequency[med_cnt] = cur_frequency[j];
                med_cnt++;
            }
        }

        return new Object[]{full_decoded_data, full_medoid_list, full_frequency, maxDecimalPlaces, dim, medoidIndexList, pageSize, kList};
    }

    public static Object[] decodeNoHuff(String filePath) throws IOException {
        final double[] POW_10 = new double[16]; 
        for (int i = 0; i < POW_10.length; i++) {
            POW_10[i] = Math.pow(10, i);
        }


        BitStreamReader reader = new BitStreamReader(filePath);

        
        int dim = (int) reader.readBits(8);
        int packSize = (int) reader.readBits(16);
        int blockSize = (int) reader.readBits(16);
        int pageCount = (int) reader.readBits(32);


        List<Integer> pageKList = new ArrayList<>();
        List<Integer> pageDataPointsList = new ArrayList<>();
        List<int[]> maxDecimalList = new ArrayList<>();
        List<long[]> minValueList = new ArrayList<>();
        List<long[]> minBaseList = new ArrayList<>();
        List<int[]> maxBaseBitLenList = new ArrayList<>();
        List<int[][]> packResMetadataList = new ArrayList<>();

        for (int i = 0; i < pageCount; i++) {
            pageKList.add((int) reader.readBits(16));
            pageDataPointsList.add((int) reader.readBits(16));

            
            int[] maxDecimals = new int[dim];
            long[] minValues = new long[dim];
            for (int j = 0; j < dim; j++) {
                maxDecimals[j] = (int) reader.readBits(8);
                int valBitLen = (int) reader.readBits(8);
                long sign = reader.readBits(1);
                long absValue = reader.readBits(valBitLen);
                minValues[j] = (sign == 0) ? absValue : -absValue;
            }
            maxDecimalList.add(maxDecimals);
            minValueList.add(minValues);

            
            long[] minBases = new long[dim];
            int[] maxBaseBits = new int[dim];
            for (int j = 0; j < dim; j++) {
                
                int valBitLen = (int) reader.readBits(8);
                long sign = reader.readBits(1);
                long absValue = reader.readBits(valBitLen);
                minBases[j] = (sign == 0) ? absValue : -absValue;
                
                maxBaseBits[j] = (int) reader.readBits(8);
            }
            minBaseList.add(minBases);
            maxBaseBitLenList.add(maxBaseBits);

            
            int packNum = (int) reader.readBits(16);
            int[][] packMetadata = new int[packNum][dim];
            for (int packIdx = 0; packIdx < packNum; packIdx++) {
                for (int j = 0; j < dim; j++) {
                    packMetadata[packIdx][j] = (int) reader.readBits(8);
                }
            }
            packResMetadataList.add(packMetadata);
        }


        

        List<long[][]> decodedMedoidsList = new ArrayList<>();
        for (int i = 0; i < pageCount; i++) {
            int k = pageKList.get(i);
            long[] minBases = minBaseList.get(i);
            int[] maxBaseBits = maxBaseBitLenList.get(i);

            long[][] medoids = new long[k][dim];
            for (int medoidIdx = 0; medoidIdx < k; medoidIdx++) {
                for (int j = 0; j < dim; j++) {
                    long offsetValue = reader.readBits(maxBaseBits[j]);
                    medoids[medoidIdx][j] = offsetValue + minBases[j];
                }
            }
            decodedMedoidsList.add(medoids);
            
        }







        

        List<long[]> decodedFrequencyList = new ArrayList<>();
        for (int i = 0; i < pageCount; i++) {
            int k = pageKList.get(i);
            long[] frequency = new long[k];
            if (k > 0) {
                int blockNum = (k + blockSize - 1) / blockSize;
                int[] blockMeta = new int[blockNum];
                for (int blockIdx = 0; blockIdx < blockNum; blockIdx++) {
                    blockMeta[blockIdx] = (int) reader.readBits(8);
                }

                long[] clusterCounts = new long[k];
                for (int blockIdx = 0; blockIdx < blockNum; blockIdx++) {
                    int start = blockIdx * blockSize;
                    int end = Math.min(start + blockSize, k);
                    int bitsToRead = blockMeta[blockIdx];
                    for (int j = start; j < end; j++) {
                        clusterCounts[j] = reader.readBits(bitsToRead);
                    }
                }
                frequency[0] = clusterCounts[0];
                for (int j = 1; j < k; j++) {
                    frequency[j] = frequency[j - 1] + clusterCounts[j];
                }
            }

            decodedFrequencyList.add(frequency);

        }


        
        List<long[][]> decodedResidualSeriesList = new ArrayList<>(); 
        for (int i = 0; i < pageCount; i++) {
            int pageDataPoints = pageDataPointsList.get(i);
            int[][] packMetadata = packResMetadataList.get(i); 

            long[][] residualSeries = new long[pageDataPoints][dim];
            int dataCounter = 0;
            if (pageDataPoints > 0) { 
                for (int packIdx = 0; packIdx < packMetadata.length; packIdx++) {
                    int[] bitsPerDim = packMetadata[packIdx];
                    int pointsInPack = (packIdx == packMetadata.length - 1)
                            ? (pageDataPoints - (packIdx * packSize))
                            : packSize;
                    for (int p = 0; p < pointsInPack; p++) {
                        if (dataCounter >= pageDataPoints) { 
                            System.err.println("Error: dataCounter exceeds pageDataPoints during residual decoding. Page: " + i);
                            break;
                        }
                        
                        for (int j = 0; j < dim; j++) {
                            long zigzagValue = reader.readBits(bitsPerDim[j]);

                            residualSeries[dataCounter][j] = (zigzagValue >> 1) ^ (-(zigzagValue & 1));



                        }
                        dataCounter++;
                    }
                    if (dataCounter >= pageDataPoints && packIdx < packMetadata.length -1) {
                        System.err.println("Warning: Residuals fully decoded before all packs processed. Page: " + i + ", Pack: " + packIdx);
                        break;
                    }
                }
            }
            decodedResidualSeriesList.add(residualSeries);
            
        }



        

        List<double[][]> decodedPagesAsDouble = new ArrayList<>();
        List<double[][]> finalMedoidsAsDoubleList = new ArrayList<>(); 
        int totalDataPoints = 0;
        int totalMedoids = 0;

        for (int i = 0; i < pageCount; i++) {
            
            int k = pageKList.get(i);
            int pageDataPoints = pageDataPointsList.get(i);
            long[][] medoids = decodedMedoidsList.get(i);

            long[] frequency = decodedFrequencyList.get(i);


            long[][] residualSeries = decodedResidualSeriesList.get(i); 
            long[] minValues = minValueList.get(i); 
            int[] maxDecimals = maxDecimalList.get(i); 

            int flength = frequency.length;
            long[] cumFrequency = new long[flength];
            cumFrequency[0]  = frequency[0];
            for(int f=1; f<flength; f++){
                cumFrequency[f] = cumFrequency[f-1] + frequency[f];
            }



            double[][] finalPageData = new double[pageDataPoints][dim];
            int currentMedoidIndex = 0;

            for (int dataIdx = 0; dataIdx < pageDataPoints; dataIdx++) {
                if (k > 0 && currentMedoidIndex < k - 1 && dataIdx >= cumFrequency[currentMedoidIndex]) {
                    currentMedoidIndex++;
                }
                long[] basePoint = (k > 0 && medoids.length > 0 && currentMedoidIndex < medoids.length) ? medoids[currentMedoidIndex] : new long[dim];

                for (int j = 0; j < dim; j++) {
                    long integerValue = basePoint[j] + residualSeries[dataIdx][j] + minValues[j];
                    finalPageData[dataIdx][j] = (double) integerValue / POW_10[maxDecimals[j]]; 
                }
            }
            decodedPagesAsDouble.add(finalPageData);
            totalDataPoints += pageDataPoints;

            double[][] pageDoubleMedoids = new double[k][dim];
            for (int medoidIdx = 0; medoidIdx < k; medoidIdx++) {
                for (int j = 0; j < dim; j++) {
                    long integerValue = medoids[medoidIdx][j] + minValues[j];
                    pageDoubleMedoids[medoidIdx][j] = (double) integerValue / POW_10[maxDecimals[j]]; 
                }
            }
            finalMedoidsAsDoubleList.add(pageDoubleMedoids);
            totalMedoids += k;
        }




        
        double[][] finalDecodedData = new double[totalDataPoints][dim];
        int currentRow = 0;
        for (double[][] pageData : decodedPagesAsDouble) {
            for (int i = 0; i < pageData.length; i++) {
                
                for (int j = 0; j < dim; j++) {
                    finalDecodedData[currentRow][j] = pageData[i][j];
                }
                currentRow++;
            }
        }


        
        double[][] finalMedoids = new double[totalMedoids][dim];
        int medoidCursor = 0;
        for (double[][] pageMedoids : finalMedoidsAsDoubleList) {
            for (int i = 0; i < pageMedoids.length; i++) {
                
                for (int j = 0; j < dim; j++) {
                    finalMedoids[medoidCursor][j] = pageMedoids[i][j];
                }
                medoidCursor++;
            }
        }

        
        long[] finalFrequencies = new long[totalMedoids];
        medoidCursor = 0; 
        for (int i = 0; i < pageCount; i++) {
            long[] pageFrequencies = decodedFrequencyList.get(i);
            for (int j = 0; j < pageFrequencies.length; j++) {
                
                finalFrequencies[medoidCursor + j] = pageFrequencies[j];
            }
            medoidCursor += pageKList.get(i);
        }




        return new Object[]{finalDecodedData, finalMedoidsAsDoubleList, decodedFrequencyList, finalMedoids, finalFrequencies, };
    }



    public static Object[] decode_PageByPage(String filePath) throws IOException {

        final double[] POW_10 = new double[16];
        for (int i = 0; i < POW_10.length; i++) {
            POW_10[i] = Math.pow(10, i);
        }

        BitStreamReader reader = new BitStreamReader(filePath);

        // =====================================================================
        // 1. 读取全局文件头 (无问题)
        // =====================================================================
        int dim = (int) reader.readBits(8);
        int packSize = (int) reader.readBits(16);
        int blockSize = (int) reader.readBits(16);
        int pageCount = (int) reader.readBits(32);
        int pageSize = (int) reader.readBits(16);

        // =====================================================================
        // 2. 预分配与初始化 (无问题)
        // =====================================================================
        long maxTotalDataPoints = (long) pageCount * pageSize;
        double[][] finalDecodedData = new double[(int) maxTotalDataPoints][dim];
        int dataCursor = 0;

        List<double[][]> finalMedoidsAsDoubleList = new ArrayList<>(pageCount);
        List<long[]> decodedFrequencyList = new ArrayList<>(pageCount);

        // =====================================================================
        // 3. 进入主解码循环
        // =====================================================================
        for (int i = 0; i < pageCount; i++) {

            // --- 3.1 & 3.2 & 3.3: 元数据、基准点、频率解码 ---
            int k = (int) reader.readBits(16);
            int pageDataPoints = (int) reader.readBits(16);
            int[] maxDecimals = new int[dim];
            long[] minValues = new long[dim];
            for (int j = 0; j < dim; j++) {
                maxDecimals[j] = (int) reader.readBits(8);
                int valBitLen = (int) reader.readBits(8);
                int sign = (int) reader.readBits(1);
                long absValue = reader.readBits(valBitLen);
                minValues[j] = (sign == 0) ? absValue : -absValue;
            }

            long[] minBases = new long[dim];
            int[] maxBaseBits = new int[dim];
            for (int j = 0; j < dim; j++) {
                int valBitLen = (int) reader.readBits(8);
                int sign = (int) reader.readBits(1);
                long absValue = (int) reader.readBits(valBitLen);
                minBases[j] = (sign == 0) ? absValue : -absValue;
                maxBaseBits[j] = (int) reader.readBits(8);
            }

            int packNum = (int) reader.readBits(16);
            int[][] packMetadata = new int[packNum][dim];
            for (int packIdx = 0; packIdx < packNum; packIdx++) {
                for (int j = 0; j < dim; j++) {
                    packMetadata[packIdx][j] = (int) reader.readBits(8);
                }
            }

            // --- 后续解码逻辑 (保持不变，因为前面的修正会使其恢复正常) ---
            double[][] medoidsAsDouble = new double[k][dim];
            long[][] medoidsAsLong = new long[k][dim];
            if (k > 0) {
                for (int medoidIdx = 0; medoidIdx < k; medoidIdx++) {
                    for (int j = 0; j < dim; j++) {
                        long longMedoid = minBases[j] + reader.readBits(maxBaseBits[j]);
                        medoidsAsLong[medoidIdx][j] = longMedoid;
                        long integerValue = longMedoid + minValues[j];
                        medoidsAsDouble[medoidIdx][j] = (double) integerValue / POW_10[maxDecimals[j]];
                    }
                }
            }
            System.out.println(Arrays.deepToString(medoidsAsDouble));
            System.out.println("");
            finalMedoidsAsDoubleList.add(medoidsAsDouble);

            long[] frequencies = new long[k];
            if (k > 0) {
                int blockNumFreq = (k + blockSize - 1) / blockSize;
                int[] blockMetaFreq = new int[blockNumFreq];
                for (int blockIdx = 0; blockIdx < blockNumFreq; blockIdx++) {
                    blockMetaFreq[blockIdx] = (int) reader.readBits(8);
                }
                for (int blockIdx = 0; blockIdx < blockNumFreq; blockIdx++) {
                    int start = blockIdx * blockSize;
                    int end = Math.min(start + blockSize, k);
                    for (int j = start; j < end; j++) {
                        frequencies[j] = reader.readBits(blockMetaFreq[blockIdx]);
                    }
                }
                for (int c = 1; c < k; c++) {
                    frequencies[c] = frequencies[c] + frequencies[c-1];
                }
            }
            System.out.println(Arrays.toString(frequencies));
            decodedFrequencyList.add(frequencies);

            // --- 3.4 & 3.5: 残差解码与数据重构 (无变化) ---
            long[] cumFrequency = new long[k];
            if (k > 0) {
                cumFrequency[0] = frequencies[0];
                for (int f = 1; f < k; f++) {
                    cumFrequency[f] = cumFrequency[f-1] + frequencies[f];
                }
            }
            if (pageDataPoints > 0) {
                long[][] residualSeries = new long[pageDataPoints][dim];
                int residualDataCounter = 0;
                for (int packIdx = 0; packIdx < packNum; packIdx++) {
                    int[] bitsPerDim = packMetadata[packIdx];
                    int startPoint = packIdx * packSize;
                    int endPoint = Math.min(startPoint + packSize, pageDataPoints);
                    for (int p = startPoint; p < endPoint; p++) {
                        if (residualDataCounter >= pageDataPoints) break;
                        for (int j = 0; j < dim; j++) {
                            long zigzagValue = reader.readBits(bitsPerDim[j]);
                            residualSeries[residualDataCounter][j] = (zigzagValue >>> 1) ^ -(zigzagValue & 1);
                        }
                        residualDataCounter++;
                    }
                }

                int currentMedoidIndex = 0;
                for (int dataIdx = 0; dataIdx < pageDataPoints; dataIdx++) {
                    if (k > 0 && currentMedoidIndex < k - 1 && dataIdx >= cumFrequency[currentMedoidIndex]) {
                        currentMedoidIndex++;
                    }
                    long[] basePoint = (k > 0) ? medoidsAsLong[currentMedoidIndex] : new long[dim];
                    for (int j = 0; j < dim; j++) {
//                        long integerValue = basePoint[j] + residualSeries[dataIdx][j];
                        long integerValue = basePoint[j] + residualSeries[dataIdx][j] + minValues[j];
                        finalDecodedData[dataCursor][j] = (double) integerValue / POW_10[maxDecimals[j]];
                    }
                    dataCursor++;
                }
            }
        }
        reader.close();

        // =====================================================================
        // 4. 最终数据修剪 (无问题)
        // =====================================================================
        double[][] trulyFinalDecodedData;
        if (dataCursor < maxTotalDataPoints) {
            trulyFinalDecodedData = new double[dataCursor][dim];
            System.arraycopy(finalDecodedData, 0, trulyFinalDecodedData, 0, dataCursor);
        } else {
            trulyFinalDecodedData = finalDecodedData;
        }

        return new Object[]{trulyFinalDecodedData, finalMedoidsAsDoubleList, decodedFrequencyList};
    }

    public static Object[] decode_PageByPage_SIMD(String filePath) throws IOException {
        
        
        
        final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;
        final int VECTOR_LANE_COUNT = SPECIES.length();

        final double[] POW_10 = new double[19];
        for (int i = 0; i < POW_10.length; i++) {
            POW_10[i] = Math.pow(10, i);
        }

        BitStreamReader reader = new BitStreamReader(filePath);

        
        
        
        int dim = (int) reader.readBits(8);
        int packSize = (int) reader.readBits(16);
        int blockSize = (int) reader.readBits(16);
        int pageCount = (int) reader.readBits(32);
        int pageSize = (int) reader.readBits(16);
        
        long maxTotalDataPoints = (long) pageCount * pageSize;
        double[][] finalDecodedData = new double[(int) maxTotalDataPoints][dim];
        int dataCursor = 0;

        List<double[][]> finalMedoidsAsDoubleList = new ArrayList<>(pageCount);
        List<long[]> decodedFrequencyList = new ArrayList<>(pageCount);

        double[] tempDivisors = new double[VECTOR_LANE_COUNT];

        
        for (int i = 0; i < pageCount; i++) {

            
            
            int k = (int) reader.readBits(16);
            int pageDataPoints = (int) reader.readBits(16);
            int[] maxDecimals = new int[dim];
            long[] minValues = new long[dim];
            for (int j = 0; j < dim; j++) {maxDecimals[j] = (int) reader.readBits(8);int valBitLen = (int) reader.readBits(8);minValues[j] = (reader.readBits(1) == 0) ? reader.readBits(valBitLen) : -reader.readBits(valBitLen);}
            long[] minBases = new long[dim];
            int[] maxBaseBits = new int[dim];
            for (int j = 0; j < dim; j++) {int valBitLen = (int) reader.readBits(8);minBases[j] = (reader.readBits(1) == 0) ? reader.readBits(valBitLen) : -reader.readBits(valBitLen);maxBaseBits[j] = (int) reader.readBits(8);}
            int packNum = (int) reader.readBits(16);
            int[][] packMetadata = new int[packNum][dim];
            for (int packIdx = 0; packIdx < packNum; packIdx++) {for (int j = 0; j < dim; j++) {packMetadata[packIdx][j] = (int) reader.readBits(8);}}
            double[][] medoidsAsDouble = new double[k][dim];
            long[][] medoidsAsLong = new long[k][dim];
            if (k > 0) {for (int medoidIdx = 0; medoidIdx < k; medoidIdx++) {for (int j = 0; j < dim; j++) {long longMedoid = minBases[j] + reader.readBits(maxBaseBits[j]);medoidsAsLong[medoidIdx][j] = longMedoid;long integerValue = longMedoid + minValues[j];medoidsAsDouble[medoidIdx][j] = (double) integerValue / POW_10[maxDecimals[j]];}}}
            finalMedoidsAsDoubleList.add(medoidsAsDouble);
            long[] frequencies = new long[k];
            if (k > 0) {int blockNumFreq = (k + blockSize - 1) / blockSize;int[] blockMetaFreq = new int[blockNumFreq];for (int blockIdx = 0; blockIdx < blockNumFreq; blockIdx++) {blockMetaFreq[blockIdx] = (int) reader.readBits(8);}for (int blockIdx = 0; blockIdx < blockNumFreq; blockIdx++) {int start = blockIdx * blockSize;int end = Math.min(start + blockSize, k);for (int j = start; j < end; j++) {frequencies[j] = reader.readBits(blockMetaFreq[blockIdx]);}}}
            decodedFrequencyList.add(frequencies);

            
            long[] cumFrequency = new long[k];
            if (k > 0) {cumFrequency[0] = frequencies[0];for (int f = 1; f < k; f++) {cumFrequency[f] = cumFrequency[f-1] + frequencies[f];}}

            if (pageDataPoints > 0) {
                long[][] residualSeries = new long[pageDataPoints][dim];
                int residualDataCounter = 0;
                for (int packIdx = 0; packIdx < packNum; packIdx++) {int[] bitsPerDim = packMetadata[packIdx];int startPoint = packIdx * packSize;int endPoint = Math.min(startPoint + packSize, pageDataPoints);for (int p = startPoint; p < endPoint; p++) {if (residualDataCounter >= pageDataPoints) break;for (int j = 0; j < dim; j++) {long zigzagValue = reader.readBits(bitsPerDim[j]);residualSeries[residualDataCounter][j] = (zigzagValue >>> 1) ^ -(zigzagValue & 1);}residualDataCounter++;}}

                int currentMedoidIndex = 0;
                for (int dataIdx = 0; dataIdx < pageDataPoints; dataIdx++) {
                    if (k > 0 && currentMedoidIndex < k - 1 && dataIdx >= cumFrequency[currentMedoidIndex]) {
                        currentMedoidIndex++;
                    }
                    long[] basePoint = (k > 0) ? medoidsAsLong[currentMedoidIndex] : new long[dim];

                    int loopBound = SPECIES.loopBound(dim);

                    int j = 0;
                    
                    for (; j < loopBound; j += VECTOR_LANE_COUNT) {
                        
                        LongVector baseVec = LongVector.fromArray(LongVector.SPECIES_64, basePoint, j);
                        LongVector residualVec = LongVector.fromArray(LongVector.SPECIES_64, residualSeries[dataIdx], j);
                        LongVector minValVec = LongVector.fromArray(LongVector.SPECIES_64, minValues, j);

                        for (int lane = 0; lane < VECTOR_LANE_COUNT; lane++) {
                            tempDivisors[lane] = POW_10[maxDecimals[j + lane]];
                        }
                        DoubleVector divisorVec = DoubleVector.fromArray(SPECIES, tempDivisors, 0);

                        LongVector sumLongVec = baseVec.add(residualVec).add(minValVec);
                        DoubleVector sumDoubleVec = (DoubleVector) sumLongVec.convert(VectorOperators.L2D, 0);
                        DoubleVector resultVec = sumDoubleVec.div(divisorVec);
                        
                        resultVec.intoArray(finalDecodedData[dataCursor], j);
                    }
                    
                    for (; j < dim; j++) {
                        long integerValue = basePoint[j] + residualSeries[dataIdx][j] + minValues[j];
                        finalDecodedData[dataCursor][j] = (double) integerValue / POW_10[maxDecimals[j]];
                    }

                    dataCursor++;
                }
            }
        }
        reader.close();

        double[][] trulyFinalDecodedData;
        if (dataCursor < maxTotalDataPoints) {
            trulyFinalDecodedData = new double[dataCursor][dim];
            System.arraycopy(finalDecodedData, 0, trulyFinalDecodedData, 0, dataCursor);
        } else {
            trulyFinalDecodedData = finalDecodedData;
        }

        return new Object[]{trulyFinalDecodedData, finalMedoidsAsDoubleList, decodedFrequencyList};
    }


    public static void main(String[] args) throws IOException {

        // adjust dataset_indx for testing
        // have to compress first, then try to decompress
        // we prepared the compressed file from ACluster in the resources, just for quick tests
        int dataset_indx = 0;
        String[] datasets= {"SSD-bench" ,"profile-income" ,"Blockchain-tr", "Crop","gas", "TY_display"};
        String dataset = datasets[dataset_indx];

        int page[] = {1000,2000,5000,10000,20000};
        int[] data_size = {8926,14825,99999, 24000,13910};

        long start_time = System.nanoTime();
        Object[] results = decode_PageByPage("src/test/resources/" + dataset + "_compressed.csv");
        long end_time = System.nanoTime();
        System.out.println("Decoding time: " + (end_time-start_time)/data_size[dataset_indx]);
    }
}
