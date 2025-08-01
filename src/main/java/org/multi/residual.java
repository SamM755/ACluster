package org.multi;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class residual {
    private long[][] residualSeries;
    private long[] globalRMin;

    
    public residual(long[][] residualSeries, long[] globalRMin) {
        this.residualSeries = residualSeries;
        this.globalRMin = globalRMin;
    }

    public static long[][] residualCalculationZigzagNoHuff_sorted(long[][] data, long[][] medoids, int[] clusterAssignment, long[] frequency,
                                                           String residualOutputFile, boolean save, int dim) {
        int n = data.length;
        int k = medoids.length;

        if (n == 0) {
            return new long[0][dim];
        }

        
        int[] writePointers = new int[k];
        int cumulativeCount = 0;
        for (int i = 0; i < k; i++) {
            writePointers[i] = cumulativeCount;
            cumulativeCount += (int) frequency[i];
        }

        long[][] sortedResidualSeries = new long[n][dim];
        for (int i = 0; i < n; i++) {
            int clusterId = clusterAssignment[i];
            long[] base = medoids[clusterId];

            
            int targetIndex = writePointers[clusterId];

            
            for (int j = 0; j < dim; j++) {
                long residual = data[i][j] - base[j];
                sortedResidualSeries[targetIndex][j] = zigzagEncode(residual);
            }

            writePointers[clusterId]++;
        }

        
        if (save) {
            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(residualOutputFile))) {
                for (long[] residualRow : sortedResidualSeries) {
                    for (int j = 0; j < dim; j++) {
                        writer.write(Long.toString(residualRow[j]));
                        if (j < dim - 1) {
                            writer.write(",");
                        }
                    }
                    writer.newLine();
                }
                writer.flush();
            } catch (IOException e) {
                System.err.println(e.getMessage());
                
            }
        }

        return sortedResidualSeries;
    }

    
    public long[][] getResidualSeries() {
        return residualSeries;
    }

    
    public long[] getGlobalRMin() {
        return globalRMin;
    }

    
    public static residual residualCalculation(long[][] data, long[][] medoids, int[] clusterAssignment, int dim) {
        long[][] residualSeries = new long[data.length][dim];
        long[] globalRMin = new long[dim];
        for(int i=0; i<dim; i++){
            globalRMin[i] = Long.MAX_VALUE;  
        }

        
        for (int valueId = 0; valueId < data.length; valueId++) {
            long[] base = medoids[clusterAssignment[valueId]];
            for(int i=0; i<dim; i++){
                residualSeries[valueId][i] = data[valueId][i] - base[i];
                if (residualSeries[valueId][i] < globalRMin[i]) {
                    globalRMin[i] = residualSeries[valueId][i];
                }
            }

            

        }

        
        for (int i = 0; i < residualSeries.length; i++) {
            for(int j=0; j<dim; j++){
                residualSeries[i][j] -= globalRMin[j];
            }
        }

        
        return new residual(residualSeries, globalRMin);
    }

    
    public static long[][] residualCalculationZigzag(long[][] data, long[][] medoids, int[] clusterAssignment, String residualOutputFile, boolean save, int dim) throws IOException {
        long[][] residualSeries = new long[data.length][dim];

        


                for (int valueId = 0; valueId < data.length; valueId++) {
                    long[] base = medoids[clusterAssignment[valueId]];
                    for(int i=0; i<dim; i++){


                        long res = data[valueId][i] - base[i];

                        residualSeries[valueId][i] = zigzagEncode(res);
                    }

                }

        
        if (save) {
            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(residualOutputFile))) {
                for(int i=0; i<residualSeries.length; i++){
                    for(int j=0; j<dim; j++){
                        writer.write(medoids[clusterAssignment[i]]+",");
                    }
                    for(int j=0; j<dim; j++){
                        writer.write(residualSeries[i][j] + ", ");
                    }
                    writer.write("\n");  
                }
                writer.flush();  
            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
        }

        return residualSeries;
    }

    
    public static long zigzagEncode(long n) {
        return (n << 1) ^ (n >> 63);  
    }

    
    public static long zigzagDecode(long encoded) {
        return (encoded >> 1) ^ (-(encoded & 1));  
    }

    public static void main(String[] args) throws IOException {
        
        long[][] data = {
                {10, 20},
                {30, 40},
                {50, 60},
                {70, 80}
        };
        long[][] medoids = {
                {10, 20},
                {35, 45},
                {55, 65}
        };
        int[] clusterAssignment = {0, 1, 2, 2};

        
        long[][] zigzagResiduals = residualCalculationZigzag(data, medoids, clusterAssignment, "", false, 2);
        System.out.println("Zigzag Residuals: " + Arrays.deepToString(zigzagResiduals));
    }
}
