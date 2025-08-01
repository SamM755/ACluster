package org.multi;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;

public class Preprocessor {

    public static int[] findMaxDecimalPlaces(double[][] data, int dim) {
        
        int[] maxDecimalPlaces = new int[dim];
        if (data == null || data.length == 0) {
            return maxDecimalPlaces;
        }

        
        for (int row = 0; row < data.length; row++) {
            
            for (int col = 0; col < dim; col++) {
                String numberStr = Double.toString(data[row][col]);

                
                if (numberStr.contains(".")) {
                    
                    String decimalPart = numberStr.split("\\.")[1];

                    
                    
                    int currentDecimalPlaces = decimalPart.equals("0") ? 0 : decimalPart.length();

                    
                    maxDecimalPlaces[col] = Math.max(maxDecimalPlaces[col], currentDecimalPlaces);
                }
                
            }
        }
        return maxDecimalPlaces;
    }
    public static int[] findMaxDecimalPlacesFromStrings(String[] data, int dim) {
        int[] maxDecimalPlaces = new int[dim];

        for (String rowStr : data) {
            String[] values = rowStr.split(";");
            if (values.length != dim) {
                throw new IllegalArgumentException("Dimension mismatch: expected " + dim + ", got " + values.length);
            }

            for (int col = 0; col < dim; col++) {
                
                String val = values[col].trim();
                
                int dotIndex = val.indexOf('.');
                if (dotIndex >= 0) {
                    
                    int decimalLen = val.length() - (dotIndex + 1);
                    maxDecimalPlaces[col] = Math.max(maxDecimalPlaces[col], decimalLen);
                }
            }
        }

        return maxDecimalPlaces;
    }

    
    public static long[][] convertToIntegers(double[][] data, int[] maxDecimalPlaces) {
        if (data == null || data.length == 0 || maxDecimalPlaces == null) {
            throw new IllegalArgumentException("Input data or maxDecimalPlaces cannot be null or empty.");
        }

        
        int numColumns = maxDecimalPlaces.length;
        for (double[] row : data) {
            if (row.length != numColumns) {
                throw new IllegalArgumentException("Row length does not match the number of dimensions.");
            }
        }

        
        long[][] integerData = new long[data.length][numColumns];

        
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < numColumns; j++) {

                integerData[i][j] = Math.round(data[i][j] * Math.pow(10, maxDecimalPlaces[j]));
            }
        }

        return integerData;
    }

    public static long[] preprocessDataWithDelta(long[] data) {
        long[] deltaData = new long[data.length];
        deltaData[0] = data[0];
        for (int i = 1; i < data.length; i++) {
            long deltaValue = data[i] - data[i - 1];
            deltaData[i] = deltaValue;
        }
        return deltaData;
    }

    
    public static void printStatsComparison(double deltaValue, double originalValue, String message) {
        if (deltaValue >= originalValue) {
            System.out.println(message);
        }
    }

    
    public static double getMedianBitWidth(long[] data) {
        return Arrays.stream(data)
                .mapToDouble(value -> Math.ceil(Math.log(Math.abs(value) + 1) / Math.log(2)))
                .sorted()
                .skip(data.length / 2)
                .findFirst()
                .orElse(0);
    }

    
    public static double getMedianAbs(long[] data) {
        return Arrays.stream(data)
                .mapToDouble(Math::abs)
                .sorted()
                .skip(data.length / 2)
                .findFirst()
                .orElse(0);
    }

    
    public static double getIQR(long[] data) {
        double q1 = getPercentile(data, 25);
        double q3 = getPercentile(data, 75);
        return q3 - q1;
    }

    
    public static double getVariance(long[] data) {
        double mean = Arrays.stream(data).average().orElse(0.0);
        return Arrays.stream(data).mapToDouble(d -> Math.pow(d - mean, 2)).average().orElse(0.0);
    }

    
    public static double getEntropy(long[] data) {
        Map<Long, Long> freqMap = new HashMap<>();
        for (long value : data) {
            freqMap.put(value, freqMap.getOrDefault(value, 0L) + 1);
        }
        double entropy = 0;
        long size = data.length;
        for (long count : freqMap.values()) {
            double p = (double) count / size;
            entropy -= p * Math.log(p) / Math.log(2);
        }
        return entropy;
    }

    
    public static double getPercentile(long[] data, int percentile) {
        int index = (int) Math.ceil(percentile / 100.0 * data.length);
        return Arrays.stream(data).sorted().skip(index).findFirst().orElse(0);
    }

    public static long[][] processDelta(long[][] data, int dim, ArrayList<Integer> dTime) {
        long[][] deltaData = new long[data.length][dim];
        for(int j=0; j<dim; j++){
            deltaData[0][j] = data[0][j];
        }
        for(int i=1; i<data.length; i++){
            long[] deltaValue = new long[dim];
            for(int j=0; j<dim; j++){
                deltaValue[j] = data[i][j] - data[i-1][j];
                deltaData[i][j] = deltaValue[j];
            }
        }
        return deltaData;
    }

    public static boolean shouldApplyDeltaForCompression(long[] data) {
        
        long[] deltaData = preprocessDataWithDelta(data);

        
        double originalMedianBitWidth = getMedianBitWidth(data);
        double deltaMedianBitWidth = getMedianBitWidth(deltaData);

        double originalAbsMedian = getMedianAbs(data);
        double deltaAbsMedian = getMedianAbs(deltaData);

        double originalIQR = getIQR(data);
        double deltaIQR = getIQR(deltaData);

        double originalVariance = getVariance(data);
        double deltaVariance = getVariance(deltaData);

        double originalEntropy = getEntropy(data);
        double deltaEntropy = getEntropy(deltaData);

        
        if (deltaMedianBitWidth < originalMedianBitWidth &&
                deltaAbsMedian < originalAbsMedian &&
                deltaIQR <= originalIQR &&
                deltaVariance <= originalVariance &&
                deltaEntropy <= originalEntropy) {
            System.out.println("适合应用delta编码");
            return true;
        } else {
            
            if (deltaMedianBitWidth >= originalMedianBitWidth) {
                System.out.println("Delta median bit width is worse.");
            }
            if (deltaAbsMedian >= originalAbsMedian) {
                System.out.println("Delta abs median is worse.");
            }
            if (deltaIQR > originalIQR) {
                System.out.println("Delta IQR is worse.");
            }
            if (deltaVariance > originalVariance) {
                System.out.println("Delta variance is worse.");
            }
            if (deltaEntropy > originalEntropy) {
                System.out.println("Delta entropy is worse.");
            }

            System.out.println("不适合应用delta编码");
            return false;
        }
    }

    public static Object[] preprocessing(String inputFilePath, int dim) throws IOException {
        String[] stringData = support.readDataAsString(inputFilePath);


        double[][] double_data = support.convertToDouble(stringData, dim);

        int[] maxDecimalPlaces = findMaxDecimalPlacesFromStrings(stringData,dim);

        long[][] convertedData = convertToIntegers(double_data, maxDecimalPlaces);

        return new Object[]{convertedData, maxDecimalPlaces};
    }

    public static void main(String[] args) {


    }

    public static Object[] convertToIntegersAndMin(double[][] data, int[] maxDecimalPlaces) {
        if (data == null || data.length == 0 || maxDecimalPlaces == null) {
            throw new IllegalArgumentException("Input data or maxDecimalPlaces cannot be null or empty.");
        }

        int numRows = data.length;
        int numColumns = maxDecimalPlaces.length;

        
        long[][] integerData = new long[numRows][numColumns];
        long[] minValues = new long[numColumns];
        Arrays.fill(minValues, Long.MAX_VALUE);

        double[] multipliers = new double[numColumns];
        for (int j = 0; j < numColumns; j++) {
            multipliers[j] = Math.pow(10, maxDecimalPlaces[j]);
        }

        for (int i = 0; i < numRows; i++) {
            if (data[i] == null || data[i].length != numColumns) {
                throw new IllegalArgumentException("Row " + i + " has incorrect number of dimensions.");
            }
            for (int j = 0; j < numColumns; j++) {
                
                long scaledValue = Math.round(data[i][j] * multipliers[j]);
                integerData[i][j] = scaledValue;

                
                minValues[j] = Math.min(minValues[j], scaledValue);
            }
        }

        return new Object[]{integerData, minValues, numRows};
    }

    public static Object[] preprocessingPos(double[][] rawData, int dim) {

        int[] maxDecimalPlaces = findMaxDecimalPlaces(rawData,dim);

        Object[] res = convertToIntegersAndMin(rawData, maxDecimalPlaces);
        int numCol = dim;
        long[][] data = (long[][]) res[0];
        long[] minValues = (long[]) res[1];
        int numRows = (int) res[2];
        long[][] preprocessed_data = new long[numRows][numCol];
        for(int i=0; i<numRows; i++){
            for(int j=0; j<numCol; j++){
                preprocessed_data[i][j] = data[i][j] - minValues[j];
            }
        }
        return new Object[]{preprocessed_data, maxDecimalPlaces, minValues, numRows};
    }


}
